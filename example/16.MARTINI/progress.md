# Progress Log

## 2026-02-11
- Initialized required workflow files: `task_plan.md`, `findings.md`, `progress.md`.
- Prepared to reproduce failure for `run_sim_1ubq.sh`.
- Ran `bash run_sim_1ubq.sh`; initial failure was environment invocation from wrong path (`.venv/bin/activate` not found when pre-sourcing in current dir). Switched to running script directly in-place since script already sources `../../source.sh` and `../../.venv/bin/activate`.
- Reproduced workflow blocker at stage handoff:
  `set_initial_position.py` failed with `IndexError: Index (-1) out of range for empty dimension` when reading `f['/output/pos'][-1, ...]`.
- Modified `set_initial_position.py`:
  - Added guard to use `/output/pos` only when it exists and has at least one frame.
  - Added fallback to `/input/pos` for empty-output files (e.g., minimize-only files).
- Re-ran workflow with short verification run:
  `NPT_EQUIL_STEPS=20 NPT_PROD_STEPS=20 FRAME_INTERVAL=5 bash run_sim_1ubq.sh`.
- Verification results:
  - `outputs/martini_test/logs/npt_equilibration.log` shows completion with checkpoint write.
  - `outputs/martini_test/logs/npt_production.log` shows simulation startup (`0 / 20 ...` progress line).
  - Checkpoints produced: `1ubq.prepared.up`, `1ubq.minimized.up`, `1ubq.npt_equil.up`, `1ubq.npt_equil_reduced.up`, `1ubq.npt_prod.up`.

## 2026-02-13
- Investigated scan sweep mismatch between `scan_lipid_diffusion.py` and `scan_water_diffusion.py`.
- Confirmed lipid previously used a narrow hardcoded grid (`T=0.75..1.00`, `tau=[1,2,4,8]`) while water used a broad scan grid.
- Updated `scan_lipid_diffusion.py` with minimal parity changes:
  - Added argparse-based entrypoint style matching water scan script.
  - Changed `build_temperature_range()` to `np.arange(0.600, 1.001, 0.02)`.
  - Changed `build_tau_range()` to merged low/high/special tau logic used by water scan.
  - Kept bilayer-specific runtime script generation and diffusion-analysis logic unchanged.
- Verification run:
  - Command: `source ../../.venv/bin/activate && source ../../source.sh && python3 scan_lipid_diffusion.py --base-dir lipid_diffusion_check`
  - Result: `Generated 651 simulation scripts.`
  - Verified `lipid_diffusion_check/tasks.txt` has 651 lines and `lipid_diffusion_check/run_scan.slurm` contains `#SBATCH --array=1-651`.

## 2026-02-13 (Job 2: Protein-in-Bilayer Packing)
- Captured user requirements for new workflow:
  - Use martinize2-style coarse-graining algorithms for MARTINI 2.2 protein model.
  - Add ions consistent with 0.15 M NaCl.
  - Allow simulation box resizing after packing.
- Drafted and stored implementation plan in `task_plan.md` for review before code changes.
- Expanded hybrid architecture requirements:
  - Production-only activation of Upside<->dry-MARTINI coupling.
  - Protein rigid during pre-production.
  - MARTINI BB computed from COM of Upside `N,CA,C,O` atoms (not vice versa).
  - Intra-protein MARTINI interactions disabled entirely.
  - Rotamer mapping used only for sidechain-environment coupling.
- Replaced `task_plan.md` with detailed phased plan for schema, data prep, Upside runtime modifications, and validation.
- Completed Phase 1 planning artifact: finalized proposed `.up` hybrid schema and execution boundaries (stage-gated activation, intra-protein MARTINI exclusions, BB-from-backbone COM mapping).
- Added new requested item to plan: optional rigid-body drift mitigation for protein (frame-to-frame RMSD alignment/rotation removal mode), with explicit need to benchmark side effects and define where alignment is applied.
- Decision confirmed: rigid-body alignment mode will affect coupling coordinates only; raw integrated coordinates remain unchanged in saved outputs.
- Implemented Phase 2 data-preparation script: `prepare_hybrid_system.py`.
  - Added CLI/config (`Config` dataclass) for protein/bilayer inputs, optional martinize command, packing cutoffs, salt concentration, box padding, and seed.
  - Implemented PDB parsing/writing and protein insertion into bilayer frame.
  - Implemented clash-based lipid pruning, box resizing, ion placement (neutralization + 0.15 M NaCl target), and packed MARTINI output writing.
  - Implemented hybrid mapping export artifacts:
    - `hybrid_mapping.h5` with `/input/hybrid_control`, `/input/hybrid_bb_map`, `/input/hybrid_env_topology`, and scaffold `/input/hybrid_sc_map`.
    - `hybrid_bb_map.json` and `hybrid_prep_summary.json`.
- Fixed runtime issues during validation:
  - NumPy 2.0 incompatibility (`np.string_`) switched to byte-string attrs in HDF5.
  - Corrected protein-CG extraction to avoid treating full-system MARTINI file as protein-only input.
- Validation runs:
  - `source ../../.venv/bin/activate && source ../../source.sh && python3 prepare_hybrid_system.py --protein-pdb pdb/1rkl.pdb --protein-cg-pdb pdb/1rkl.MARTINI.pdb --bilayer-pdb pdb/bilayer.MARTINI.pdb --output-dir outputs/hybrid_1rkl_test2`
  - Output summary reports: protein CG atoms `68`, total atoms `1263`, BB map entries `31`, Na/Cl added `87/84` (protein charge `-3`, salt pairs `84`).
  - HDF5 sanity checks passed for dataset presence and shapes.
- Implemented Phase 3 (Upside hybrid input parsing + stage-gated flags) with minimal non-core source changes:
  - `src/martini.cpp`:
    - Added `martini_hybrid` runtime registry and HDF5 parser for `/input/hybrid_control`, `/input/hybrid_bb_map`, `/input/hybrid_env_topology`, and schema-check path for `/input/hybrid_sc_map`.
    - Added validation for required datasets/shapes and atom-count consistency (`protein_membership` length must equal `n_atom`).
    - Added stage-gated runtime state fields (`enabled`, `activation_stage`, `active`, `preprod_mode`, `exclude_intra_protein_martini`) and query helpers.
    - Hooked hybrid active-flag updates into `martini_stage_params::switch_simulation_stage`.
    - Updated stage-parameter default stage handling to read optional `/input/stage_parameters` attribute `current_stage` (default `production`).
  - Minimal wiring outside `martini.cpp`:
    - `src/main.h`: added forward declarations for `martini_hybrid` API.
    - `src/main.cpp`: registered hybrid config for each engine after stage-parameter registration.
- Build verification:
  - Command: `cmake --build ../../obj -j4`
  - Result: build completed successfully for `upside`, `upside_engine`, and `upside_calculation`.
- Implemented Phase 4 hybrid coupling scaffolding in Upside with minimal source footprint (core changes concentrated in `src/martini.cpp`):
  - Added hybrid runtime state sharing by engine and by coordinate node for use inside `MartiniPotential`.
  - Added BB mapping payload loading (`bb_residue_index`, `atom_indices`, `atom_mask`, `weights`, `bb_atom_index`) plus fallback BB index inference from `/input/residue_ids` + `/input/atom_names|atom_roles`.
  - Added env membership loading and active intra-protein interaction exclusion checks in MARTINI pair loop.
  - Added BB COM refresh before MARTINI pair-force evaluation.
  - Added BB-gradient projection to mapped atoms after MARTINI pair-force evaluation and zeroing of BB gradient.
  - Kept coupling stage-gated by existing stage switch (`martini_stage_params::switch_simulation_stage`).
- Updated data-prep export for Phase 4 compatibility:
  - `prepare_hybrid_system.py` now writes `/input/hybrid_bb_map/bb_atom_index`.
- Build and validation:
  - Rebuilt C++ targets successfully: `cmake --build ../../obj -j4`.
  - Regenerated hybrid mapping file and verified `bb_atom_index` dataset exists and is populated.
- Implemented Phase 5 sidechain-environment mapping path with minimal integration footprint:
  - `src/martini.cpp`:
    - Extended `hybrid_sc_map` parsing to load runtime projection payloads:
      `proxy_atom_index`, `rotamer_probability`, `proj_target_indices`, `proj_weights`.
    - Added sidechain gradient projection function that maps proxy gradients to Upside target atoms and clears proxy gradients.
    - Integrated sidechain projection into `MartiniPotential::compute_value` after MARTINI pair force accumulation and before BB projection.
    - Reused existing intra-protein MARTINI pair exclusion so sidechain proxies only receive protein-environment MARTINI interactions in hybrid-active stage.
  - `prepare_hybrid_system.py`:
    - Added `collect_sc_map()` to export deterministic sidechain proxy mapping from MARTINI sidechain beads to residue `N,CA,C,O` targets.
    - Updated HDF5 export to write non-empty `/input/hybrid_sc_map` datasets including `proxy_atom_index` and `rotamer_probability`.
    - Added `sc_map_entries` to preparation summary.
- Validation:
  - Rebuild succeeded: `cmake --build ../../obj -j4`.
  - Regenerated mapping and verified sidechain datasets and counts:
    - `rotamer_id` shape `(37,)`, `proxy_atom_index` shape `(37,)`, `rotamer_probability` sum `37.0`, `proj_target_indices` shape `(37,4)`.
- Implemented Phase 6 coupling-frame rigid-body drift mitigation + diagnostics in `src/martini.cpp`:
  - Added optional hybrid control attributes:
    - `coupling_align_enable` (0/1)
    - `coupling_align_debug` (0/1)
    - `coupling_align_interval` (int)
  - Added coupling-only alignment transform builder based on BB-frame rigid transform between consecutive steps (rotation + translation), without modifying saved raw trajectories.
  - In `MartiniPotential::compute_value`, when enabled and hybrid-active:
    - Transforms protein coordinates into coupling frame for MARTINI pair-force evaluation.
    - Rotates protein gradient contributions back to raw coordinate frame before accumulation.
  - Added periodic diagnostic printout of coupling-frame drift (`rot_deg`, `trans`) controlled by debug interval.
- Validation:
  - Rebuild succeeded after implementation/fixes: `cmake --build ../../obj -j4`.
- Added default alignment-control attrs to `prepare_hybrid_system.py` HDF5 export (`/input/hybrid_control`):
  - `coupling_align_enable=0`, `coupling_align_debug=0`, `coupling_align_interval=100`.
- Verified generated mapping contains these attrs (`outputs/hybrid_1rkl_test5/hybrid_mapping.h5`).
- Completed Phase 7 close-out tasks:
  - Added hybrid mapping validator utility: `validate_hybrid_mapping.py`.
  - Added usage/testing documentation: `HYBRID_VALIDATION.md` (build checks, schema validation, runtime smoke flow, known limits).
  - Ran schema validation successfully on `outputs/hybrid_1rkl_test5/hybrid_mapping.h5` with atom-count consistency check.
  - Executed 1-step runtime smoke on synthetic hybrid `.up` input and confirmed hybrid parser/runtime activation path without crash.
- Project state: implementation phases are complete; remaining work is extended physics benchmarking on longer production trajectories (documented as post-implementation validation).
- Re-ran mapping validation from current workspace:
  - Initial attempt used wrong flag (`--h5`) and failed with CLI usage error.
  - Correct command succeeded:
    `source ../../.venv/bin/activate && source ../../source.sh && python3 validate_hybrid_mapping.py outputs/hybrid_1rkl_test5/hybrid_mapping.h5 --n-atom 1263`
  - Output summary: `n_atom=1263`, `n_protein=68`, `n_env=1195`, `n_bb=31`, `n_sc_rotamer_rows=37`.
- Implemented missing pre-production rigid enforcement for hybrid mode in `../../src/martini.cpp`:
  - Added stage-aware dynamic rigid atom set support inside `martini_fix_rigid`.
  - Added hybrid-stage hooks to apply dynamic rigid atoms before production and release at production.
  - Added automatic pre-production rigid set construction from hybrid metadata:
    entire protein (`protein_membership>=0`) + selected lipid headgroup roles.
- Updated lipid headgroup selection policy per request:
  - Default/explicit pre-production headgroup roles now `PO4` only.
  - Runtime default changed in `../../src/martini.cpp`.
  - Mapping export now writes `preprod_lipid_headgroup_roles=PO4` in `prepare_hybrid_system.py`.
- Rebuilt C++ targets successfully after changes:
  - `cmake --build ../../obj -j4`
  - Result: success for `upside`, `upside_engine`, `upside_calculation`.

## 2026-02-14 (Hybrid Double-Check Audit)
- Re-read `task_plan.md` and audited implementation coverage across:
  - `../../src/martini.cpp` hybrid parser/runtime path (activation gating, BB refresh/projection, sidechain projection, intra-protein exclusion, coupling alignment).
  - `prepare_hybrid_system.py` mapping export path.
  - `validate_hybrid_mapping.py` schema checks.
- Ran fresh mapping generation for reproducible verification:
  - `source ../../.venv/bin/activate && source ../../source.sh && python3 prepare_hybrid_system.py --protein-pdb pdb/1rkl.pdb --protein-cg-pdb pdb/1rkl.MARTINI.pdb --bilayer-pdb pdb/bilayer.MARTINI.pdb --output-dir outputs/hybrid_audit`
- Ran current validator on generated mapping:
  - `source ../../.venv/bin/activate && source ../../source.sh && python3 validate_hybrid_mapping.py outputs/hybrid_audit/hybrid_mapping.h5 --n-atom 1263`
  - Result: passes (`n_atom=1263`, `n_protein=68`, `n_env=1195`, `n_bb=31`, `n_sc_rotamer_rows=37`).
- Ran targeted index-space sanity script on generated mapping and found blocker:
  - BB projection targets: `88/124` valid mapped indices are `>= n_protein (68)`.
  - SC projection targets: `112/148` valid mapped indices are `>= n_protein (68)`.
  - Interpretation: many projection targets point outside protein region (into environment atom range for this packed system).
- Updated `task_plan.md` Known Errors / Blockers with:
  - Mapping index-space mismatch between atomistic target indices and MARTINI packed coordinate indexing.
  - Validation gap: missing bounds/protein-membership checks for BB/SC projection target indices.

## 2026-02-14 (Hybrid Mapping Fix Implementation)
- Implemented mapping index-space fix in `prepare_hybrid_system.py`:
  - `collect_bb_map()` now emits MARTINI-space BB projection targets (`atom_indices=[bb_atom_index,-1,-1,-1]`, weights `1,0,0,0`).
  - `collect_sc_map()` now emits MARTINI-space SC projection targets routed to residue BB (`proj_target_indices=[residue_bb_index,-1,-1,-1]`, weights `1,0,0,0`).
- Strengthened mapping validation in `validate_hybrid_mapping.py`:
  - Added bounds checks for BB/SC proxy and projection target indices against `protein_membership` length.
  - Added protein-local checks requiring BB/SC proxies and all valid projection targets to satisfy `protein_membership>=0`.
- Added runtime defensive checks in `../../src/martini.cpp` hybrid parser:
  - Rejects out-of-bounds or non-protein BB/SC proxy and target indices at parse time.
- Updated projection logic in `../../src/martini.cpp`:
  - `project_bb_gradient_if_active` and `project_sc_gradient_if_active` now clear proxy gradients before redistribution, preventing self-target cancellation.
- Verification:
  - Regenerated mapping:
    `source ../../.venv/bin/activate && source ../../source.sh && python3 prepare_hybrid_system.py --protein-pdb pdb/1rkl.pdb --protein-cg-pdb pdb/1rkl.MARTINI.pdb --bilayer-pdb pdb/bilayer.MARTINI.pdb --output-dir outputs/hybrid_fix`
  - Validator pass:
    `python3 validate_hybrid_mapping.py outputs/hybrid_fix/hybrid_mapping.h5 --n-atom 1263`
    -> `n_atom=1263`, `n_protein=68`, `n_env=1195`, `n_bb=31`, `n_sc_rotamer_rows=37`.
  - Target-space sanity check:
    - BB targets outside protein: `0/31`
    - SC targets outside protein: `0/37`
  - Backward check against old bad mapping:
    `python3 validate_hybrid_mapping.py outputs/hybrid_audit/hybrid_mapping.h5 --n-atom 1263`
    now fails with `ValueError: BB target index is not protein atom ...`.
  - Rebuild succeeded after C++ changes:
    `cmake --build ../../obj -j4`
    (warnings unchanged from existing code paths).

## 2026-02-14 (Mass/Weight Audit)
- Re-read `task_plan.md` and traced particle-mass flow in code:
  - `prepare_martini.py` mass generation (`/input/mass`) and force-field lookup.
  - `../../src/martini.cpp` mass loading (`martini_masses::load_masses_for_engine`).
  - `../../src/deriv_engine.cpp` mass-aware integration updates.
- Confirmed implementation logic:
  - MARTINI masses are reduced by 12 at input generation (`mass[i] = ff_mass / 12.0`).
  - If `/input/mass` exists, integration uses those masses (`x += dt * p / mass` path).
  - If `/input/mass` does not exist, integrator uses unit-mass behavior.
- Verified existing generated MARTINI checkpoints:
  - `outputs/martini_test/checkpoints/bilayer.stage_*.up` have `/input/mass` all `6.0` (consistent with dry-MARTINI `72/12`).
- Verified mixed-type reference input:
  - `inputs/1rkl.up` mass values are `{6.0, 4.5, 3.0}` (all in reduced units, i.e., divided by 12 from source FF masses).
- Attempted protein-inclusive regeneration for direct dry-MARTINI hybrid mass validation:
  - Command: `python3 prepare_martini.py 1rkl --stage minimization outputs/masscheck_1rkl`
  - Current blocker: fails on unknown DOPC atom `D2A` in `pdb/1rkl.MARTINI.pdb`.
  - Additional mass-table mismatch identified: `pdb/1rkl_proa.itp` uses MARTINI3-style bead types (`TC*`,`TN*`,`SQ*`, etc.) not present in `ff_dry/dry_martini_v2.1.itp` mass table.

## 2026-02-14 (Probabilistic Rotamer Force Path + Single-Rotamer C++ Removal)
- Re-read `task_plan.md` and aligned implementation to remove deterministic/single-rotamer sidechain projection path in C++.
- Updated `../../src/martini.cpp`:
  - `project_sc_gradient_if_active(...)` now consumes only per-rotamer row gradients and projects them grouped by proxy (`sc_rows_by_proxy`).
  - Removed old pointer/null fallback call behavior in `compute_value`; SC projection call now always uses row-gradient vector.
  - Fixed multi-row/proxy projection correctness by clearing proxy force once per proxy before accumulating all rotamer-row projections.
- Build/verification:
  - Command: `source ../../.venv/bin/activate && source ../../source.sh && cmake --build ../../obj -j4`
  - Initial sandbox attempt failed with write-permission error for `../../obj` (expected sandbox restriction).
  - Re-ran with escalation; build succeeded for `upside`, `upside_calculation`, and `upside_engine`.
  - No compile errors from the new rotamer projection changes; only pre-existing warnings remain.

## 2026-02-14 (CPP-Only Probabilistic Projection Clarification + Runtime Update)
- Applied user clarification: keep preparation scripts unchanged; implement probabilistic sidechain behavior only in Upside C++ runtime (production-stage hybrid coupling path).
- Updated `../../src/martini.cpp` probabilistic SC pair evaluation:
  - Per-rotamer SC positions are now built from mapping projection targets (`proj_target_indices`, `proj_weights`) plus per-row `local_pos`, rather than from raw MARTINI proxy atom coordinates.
  - Added `build_sc_row_proxy_pos(...)` in `MartiniPotential::compute_value` to construct coupling-frame row positions directly from protein mapping targets.
  - Removed remaining raw-proxy dependency in the probabilistic branch (only the MARTINI pair parameters and pair membership are taken from proxy/env pair definitions).
- Rebuilt after C++ changes:
  - `source ../../.venv/bin/activate && source ../../source.sh && cmake --build ../../obj -j4`
  - Build succeeded (`upside`, `upside_engine`, `upside_calculation`), with only pre-existing warnings.

## 2026-02-14 (Sidechain Placement Correction: Upside Placement + Live Rotamer Marginals)
- Re-read `task_plan.md` and applied revised architecture in C++ runtime only (no preparation-script changes):
  - Sidechain placement should come from Upside placement node per rotamer state.
  - Probabilities should come from live Upside rotamer marginals in production.
- Updated `../../src/martini.cpp`:
  - Added hybrid runtime metadata for rotamer/placement integration:
    - detects rotamer node (`rotamer*`) and sidechain placement node (`placement*_point_vector_only` excluding `_CB`) from H5.
    - parses placement `id_seq` + `affine_residue` into per-residue/per-rotamer placement groups.
  - Added runtime row-expansion path:
    - when `hybrid_sc_map` is single-state per residue, expands rows in C++ to all available placement rotamer states for that residue.
    - keeps projection targets/weights unchanged and records row-to-placement-group mapping.
  - Added live probability refresh path:
    - reads rotamer marginals every MARTINI compute via new rotamer logs (`node_lookup`, `node_marginal`).
    - maps each SC row to `(n_rotamer,node_id,rotamer)` and updates row probability from current rotamer marginal.
    - falls back to `hybrid_sc_map/rotamer_probability` when marginals are unavailable.
  - Corrected SC coordinate source order in force loop:
    - first uses placement-group centroid from `placement*_point_vector_only` output for each rotamer row.
    - falls back to `hybrid_sc_map` projection targets + local offset when placement mapping is unavailable.
  - Added runtime scheduling dependency wiring:
    - registers rotamer/placement as parents of `martini_potential` so those nodes are computed before hybrid MARTINI coupling each step.
- Updated `../../src/rotamer.cpp`:
  - Added `get_value_by_name("node_lookup")` and `get_value_by_name("node_marginal")` for runtime access to current per-node rotamer marginals.
- Build and runtime-facing checks:
  - Build (with escalation due obj dir permissions):
    `source ../../.venv/bin/activate && source ../../source.sh && cmake --build ../../obj -j4`
    -> success (`upside`, `upside_engine`, `upside_calculation`), warnings unchanged from prior code.
  - Verified new rotamer logs through Python engine API:
    - config: `../../example/01.GettingStarted/inputs/1UBQ.up`
    - `node_lookup` shape `(76,2)`, `node_marginal` shape `(76,6)`, marginal row sums = 1.
  - Re-ran hybrid mapping validator for existing artifact:
    `python3 validate_hybrid_mapping.py outputs/hybrid_fix/hybrid_mapping.h5 --n-atom 1263`
    -> passes.
- Final compatibility adjustment:
  - Added C++ fallback mapping of existing SC rows to placement groups by `(residue_index, rotamer_id)` even when row expansion is not triggered.
  - Rebuilt again after this adjustment (`cmake --build ../../obj -j4`) and re-verified `rotamer/node_lookup` + `rotamer/node_marginal` API access via `upside_engine.py`.

## 2026-02-16 (Production MARTINI Protein-Interaction Rule Audit)
- Re-read `task_plan.md` and audited production-stage interaction gating in `../../src/martini.cpp`.
- Confirmed runtime filters are applied in all relevant paths while hybrid is active:
  - Nonbonded pair loop (`martini_potential`) via `skip_pair_if_intra_protein` and `allow_protein_pair_by_rule`.
  - Bonded pair loop (`dist_spring`) via `allow_intra_protein_pair_if_active`.
  - Multi-body bonded terms (`angle_spring`, `dihedral_spring`) via `allow_multibody_term_if_active`.
- Verified stage gating metadata in checkpoints:
  - `outputs/martini_test_1rkl_hybrid/checkpoints/1rkl.stage_6.2.up`: `current_stage=minimization`, `activation_stage=production`.
  - `outputs/martini_test_1rkl_hybrid/checkpoints/1rkl.stage_7.0.up`: `current_stage=production`, `activation_stage=production`.
- Ran direct production checkpoint census (`1rkl.stage_7.0.up`) using `/input/atom_roles`, `/input/residue_ids`, `/input/hybrid_env_topology/protein_membership`:
  - Nonbonded intra-protein pairs (`/input/potential/martini_potential/pairs`):
    - `BB-BB`: `435`, allowed by rule: `0`.
    - `SC-SC`: `552`, allowed by rule: `0`.
    - `BB-SC diff-res`: `1020`, allowed by rule: `0`.
    - `BB-SC same-res`: `6`, allowed by rule: `6`.
  - Bonded intra-protein pairs (`/input/potential/dist_spring/id`):
    - `BB-BB`: `30`, allowed by rule: `0`.
    - `SC-SC`: `9`, allowed by rule: `0`.
    - `BB-SC same-res`: `28`, allowed by rule: `28`.
  - Multi-body bonded terms:
    - `angle_spring`: intra-protein terms `0`.
    - `dihedral_spring`: intra-protein terms `3`, allowed by rule: `0` (all blocked).
- Conclusion: production-stage runtime behavior matches requested rules:
  1. No MARTINI BB-BB.
  2. No MARTINI SC-SC.
  3. MARTINI SC-BB only for same residue.

## 2026-02-17 (Production Non-Protein Hard-Sphere Enforcement)
- Re-read `task_plan.md` and audited current production-stage nonbonded behavior.
- Verified current behavior before fix:
  - Stage `7.0` had `lj_soften=0`, `coulomb_soften=0`, but still `potential_type=lj_coulomb`.
  - `1rkl.stage_7.0.up` coefficients contain regular dry-MARTINI LJ eps/sigma values; this is not hard-sphere behavior.
- Implemented runtime production-only hard-sphere-like mode in `../../src/martini.cpp`:
  - Added `hybrid_control` flag `production_nonprotein_hard_sphere` (default `1`) into hybrid runtime state parsing.
  - Added pair-mode branch in `MartiniPotential`:
    - Non-protein/non-protein pairs in hybrid-active production use repulsive-only WCA-like interaction (no Coulomb).
    - Protein-involving pairs continue using existing dry-MARTINI LJ+Coulomb path (for hybrid protein-environment coupling).
  - Added runtime parse log field `nonprotein_hs=<0|1>` for traceability.
- Updated mapping export in `prepare_hybrid_system.py`:
  - `hybrid_control` now explicitly writes `production_nonprotein_hard_sphere=1`.
- Build/test results:
  - Rebuilt successfully: `cmake --build ../../obj -j4`.
  - Production 1-step smoke (temp copy) shows hybrid parse with `nonprotein_hs=1` and stable run.
  - Control comparison with same file but `production_nonprotein_hard_sphere=0` confirms toggle is effective (`nonprotein_hs=0` in parse log and different initial potential).

## 2026-02-17 (BB Force Split by Weight Fractions)
- Implemented BB->mapped-atom force projection as normalized weight fractions in `../../src/martini.cpp`:
  - `project_bb_gradient_if_active(...)` now computes `wsum` over valid mapped targets.
  - Applies `sens[target] += (w / wsum) * sens[bb]` rather than raw `w * sens[bb]`.
  - Leaves BB gradient unchanged if no valid mapped targets exist (`wsum <= 0`) to avoid force loss.
- This matches the rigid-composite distribution law (subobject force proportional to its weight share) and is consistent with BB position reconstruction normalization.
- Rebuilt successfully after patch:
  - `source ../../.venv/bin/activate && source ../../source.sh && cmake --build ../../obj -j4`
- Sidechain-to-backbone force-transfer enforcement:
  - Confirmed current generated mapping can contain SC projection targets pointing to SC proxy indices (not BB) due key-space mismatch between AA and CG residue keys.
  - Updated `../../src/martini.cpp` to enforce runtime SC force projection onto corresponding residue MARTINI `BB` target (`sc_row_bb_target`) derived from `sc_residue_index` + `bb_residue_index/bb_atom_index`.
  - Legacy `sc_proj_target_indices/sc_proj_weights` projection is now fallback-only when residue->BB lookup is unavailable.
  - Rebuilt after enforcement (`cmake --build ../../obj -j4`) successfully.

## 2026-02-14 (Integration-Cycle Audit vs Requested Hybrid Sequence)
- Audited runtime step ordering in:
  - `../../src/main.cpp` integration loop dispatch (`integration_cycle` call sites).
  - `../../src/deriv_engine.cpp` compute/propagation and sub-step integration flow.
  - `../../src/martini.cpp` hybrid coupling internals (BB refresh, per-rotamer SC env force, SC->BB transfer, BB->backbone projection).
- Confirmed:
  - Hybrid force path matches requested SC/BB/environment transfer semantics.
  - SC probabilistic force is accumulated per rotamer row and transferred to corresponding residue BB first, then BB force is redistributed by BB map.
  - Integration uses repeated force evaluations per timestep (3 stages for Verlet/Predescu or inner-step loops for multistep), so requested sequence occurs per sub-step rather than once per full outer timestep.

## 2026-02-14 (BB Reference Metadata Export + MARTINI Scheduling Hardening)
- Re-read `task_plan.md` and implemented the pending preparation/runtime updates requested for BB reference metadata and integration ordering.
- Updated `prepare_hybrid_system.py`:
  - Added BB reference metadata export in `/input/hybrid_bb_map`:
    - `reference_atom_names` (`N,CA,C,O`)
    - `reference_atom_indices` (`n_bb x 4`)
    - `reference_atom_coords` (`n_bb x 4 x 3`)
    - `bb_comment` (per-BB string)
- Updated `validate_hybrid_mapping.py`:
  - Added optional schema checks for the new BB reference/comment datasets (shape validation only; no runtime-index enforcement on reference fields).
- Updated `../../src/martini.cpp`:
  - Added runtime dependency wiring so each `martini_potential*` node depends on all `CoordNode` outputs, ensuring MARTINI coupling evaluates after coordinate-node updates in every integration sub-step.
  - Kept explicit rotamer/placement dependencies in place.
- Build/test results:
  - Build attempt in sandbox failed due write permission in `obj/` (`Operation not permitted` on `compiler_depend.make.tmp*` files).
  - Re-ran build with escalation and succeeded:
    `source .venv/bin/activate && source source.sh && cmake --build obj -j4`
  - Generated fresh mapping artifact:
    `source ../../.venv/bin/activate && source ../../source.sh && python prepare_hybrid_system.py --protein-cg-pdb pdb/1rkl.MARTINI.pdb --output-dir outputs/hybrid_meta_cycle_check`
  - Validator pass on regenerated mapping:
    `python validate_hybrid_mapping.py outputs/hybrid_meta_cycle_check/hybrid_mapping.h5`
    -> `n_atom=1263`, `n_protein=68`, `n_env=1195`, `n_bb=31`, `n_sc_rotamer_rows=37`.

## 2026-02-15 (New End-to-End Script: run_sim_1rkl.sh)
- Added new workflow driver script: `run_sim_1rkl.sh` (modeled after `run_sim_bilayer.sh`).
- Script now automates full hybrid path from preparation to simulation to trajectory extraction:
  - Stage 0 hybrid prep via `prepare_hybrid_system.py` (packed MARTINI + `hybrid_mapping.h5` export).
  - Runtime MARTINI source file setup (`pdb/<runtime_id>.MARTINI.pdb`, matching protein ITP copy).
  - Stage-wise `.up` generation via `prepare_martini.py` for 6.0 -> 7.0.
  - Hybrid mapping injection into each stage file (`hybrid_control`, `hybrid_bb_map`, `hybrid_sc_map`, `hybrid_env_topology`).
  - Stage label patching (`minimization` pre-production, `production` at stage 7.0) so production-only hybrid coupling is enabled only in production.
  - Simulation runs and VTF extraction for stages 6.1, 6.6, and 7.0.
- Added robust path behavior:
  - Script resolves its own directory, sources env from absolute paths, then `cd`s to script directory so relative paths work regardless of caller cwd.
  - Input path checks support both relative and absolute overrides.
- Verification:
  - `bash -n run_sim_1rkl.sh` passed.
  - Set executable bit: `chmod +x run_sim_1rkl.sh`.
- Smoke run (minimal steps) status:
  - Command:
    `RUN_DIR=outputs/martini_test_1rkl_hybrid_smoke MIN_60_MAX_ITER=5 MIN_61_MAX_ITER=5 EQ_62_NSTEPS=1 EQ_63_NSTEPS=1 EQ_64_NSTEPS=1 EQ_65_NSTEPS=1 EQ_66_NSTEPS=1 PROD_70_NSTEPS=1 EQ_FRAME_STEPS=1 PROD_FRAME_STEPS=1 HYBRID_VALIDATE=1 bash run_sim_1rkl.sh`
  - Fixed one script bug found by smoke:
    - moved `set -euo pipefail` after env sourcing because `source.sh` reads unset vars under `-u`.
  - Current runtime blocker remains upstream in force-field typing:
    - `prepare_martini.py` aborts with `Mass not found for atom type 'Q5'` against current `ff_dry/dry_martini_v2.1.itp` mass table.

## 2026-02-15 (martinize2 Enforcement + Compatibility Guard, superseded)
- Updated `run_sim_1rkl.sh` to enforce MARTINI2-compatible protein topology generation for dry-MARTINI workflow:
  - Added configurable martinize controls:
    - `MARTINIZE_ENABLE` (default `1`)
    - `MARTINIZE_FF` (default `martini22`)
    - `MARTINIZE_FROM_FF` (default `charmm`)
    - `MARTINIZE_CMD_TEMPLATE` (templated command with `{input},{cg_pdb},{top},{ff},{from_ff},{molname}`)
  - Added `prepare_protein_inputs()` stage:
    - runs `martinize2.py` when enabled
    - resolves generated protein `.itp` from `.top` includes
    - validates generated/provided protein bead types against `ff_dry/dry_martini_v2.1.itp` atomtype mass table before hybrid prep proceeds.
- Fixed template-handling bug:
  - Avoided bash `${...:-...}` brace conflict by assigning `MARTINIZE_CMD_TEMPLATE` default in a separate `if` block.
- Validation:
  - `bash -n run_sim_1rkl.sh` passes.
  - With `MARTINIZE_ENABLE=0`, script now fails fast with explicit incompatible bead-type list (`Q5`, `TC*`, `TN*`, etc.) instead of failing later in `prepare_martini.py`.
  - With default `MARTINIZE_ENABLE=1`, script correctly invokes martinize2 with `-ff martini22`; current local environment still lacks `vermouth`, so execution aborts early with explicit dependency hint.
- Follow-up fixes:
  - Corrected command-template handling for martinize placeholders (`{input}`, etc.) by avoiding bash `${...:-...}` parsing conflict.
  - Relaxed upfront file validation so `PROTEIN_CG_PDB`/`PROTEIN_ITP` are required only when `MARTINIZE_ENABLE=0` (generated dynamically otherwise).

## 2026-02-15 (Switch to martinize.py for MARTINI 2.2, no vermouth)
- Replaced `run_sim_1rkl.sh` coarse-graining path from `martinize2.py` to local classic `martinize.py`:
  - `MARTINIZE_SCRIPT` default is now `${SCRIPT_DIR}/martinize.py`.
  - `run_martinize()` now runs `python3 martinize.py ... -ff martini22 -name PROA`.
  - Removed `martinize2`/`vermouth` runtime dependency from this workflow path.
- Kept dry-MARTINI compatibility guard:
  - `prepare_protein_inputs()` still resolves generated `.itp` from `.top`.
  - `assert_itp_types_have_masses()` still validates generated/provided ITP atom types against `ff_dry/dry_martini_v2.1.itp`.
- Updated `prepare_martini.py` for martinize.py output alignment:
  - Protein bead detection no longer depends only on `PROA` token in PDB lines.
  - Added residue-sequence fallback mapping so martinize ITP `resnr` can map robustly when packed-PDB residue IDs differ.
- Validation:
  - Direct command:
    `source ../../.venv/bin/activate && source ../../source.sh && python3 martinize.py -f pdb/1rkl.pdb -x outputs/martinize_py_validation/1rkl.MARTINI.pdb -o outputs/martinize_py_validation/1rkl.top -ff martini22 -name PROA`
  - Log confirms forcefield selection:
    `The martini22 forcefield will be used.`
  - Generated protein ITP (`PROA_A.itp`) contains MARTINI2-style bead types (`Qa`, `Qd`, `P*`, `C*`) and no MARTINI3-only types (`Q5`, `TC*`, `TN*`, `SQ*`).

## 2026-02-15 (Stage 6.0 HDF5 Failure Debug + Fix)
- Investigated user-reported Stage 6.0 abort:
  - `ERROR: error -1` followed by HDF5 close-ID diagnostics (`H5Gclose/H5Fclose`).
- Isolated failure scope with direct run checks:
  - `1rkl.stage_6.0.nohyb.up` (hybrid groups removed) runs successfully.
  - Any file containing `/input/hybrid_control` triggered the failure, including minimal control-only variants.
- Traced failure path into `../../src/martini.cpp` hybrid parser:
  - Failure occurred during pre-runtime node-name discovery inside `read_hybrid_settings()`.
  - Legacy HDF5 helper patterns using object-name `"."` in this path were brittle against current HDF5 behavior in this environment and produced raw `error -1` without context.
- Implemented C++ fix in `../../src/martini.cpp`:
  - Hardened `attribute_exists_hybrid()` so `"."` uses `H5Aexists(loc_id, ...)` directly (no object open/close aliasing risk).
  - Removed HDF5-time node-name discovery from `read_hybrid_settings()` (`find_rotamer_node_name` / `find_sc_placement_node_name` path).
  - Added runtime node discovery from `DerivEngine` graph in `register_hybrid_for_engine()`:
    - `discover_rotamer_node_name_from_engine()`
    - `discover_sc_placement_node_name_from_engine()` (prefers `placement_fixed_point_vector_only` / `placement_point_vector_only`, excludes `_CB`).
- Rebuilt binaries:
  - `source ../../.venv/bin/activate && source ../../source.sh && cmake --build ../../obj -j4`
  - Build succeeded.
- Verification:
  - Previously failing `1rkl.stage_6.0.prepared.up` now runs successfully through minimization.
  - Minimal full workflow smoke run now passes Stage 6.0/6.1/6.2/6.3/6.4/6.5/6.6/7.0 simulation stages with hybrid control present.
  - Remaining downstream issue is separate: `extract_martini_vtf.py` raises boolean-mask size mismatch (`1260` vs `3758`) when using non-runtime structure mapping for this run.

## 2026-02-15 (Hybrid Packing Geometry + Ion Count Corrections)
- Addressed user-reported geometry/salinity issues in `prepare_hybrid_system.py`:
  - Bilayer too small in XY for elongated protein.
  - XY box edges not matching bilayer edges.
  - Excessive NA/CL counts from full-box concentration estimate.
- Implemented packing fixes:
  - Added explicit lipid recognition helper for both `DOP` and `DOPC` (`lipid_resname`).
  - Added bilayer tiling + cropping pipeline for lipids:
    - `tile_and_crop_bilayer_lipids(...)`
    - repeats bilayer lattice using CRYST1 tile lengths (fallback: lipid extents),
    - crops by residue COM to target XY window.
  - Target XY window now enforces at least one protein-length margin beyond each side of protein:
    - `target_min = protein_min - protein_span - box_padding_xy`
    - `target_max = protein_max + protein_span + box_padding_xy`
  - Added XY box-edge lock to bilayer edges:
    - `set_box_from_lipid_xy(...)` sets `box_x/box_y` from kept lipid extents and shifts all atoms so lipid min XY is exactly `0,0`.
  - Removed global coordinate modulo wrap for all atoms (it was folding boundary lipids and breaking edge alignment).
- Implemented ion-count correction:
  - Added template calibration for effective ion-accessible volume:
    - `infer_effective_ion_volume_fraction_from_template(...)`
  - Updated salt-pair estimate to scale by this fraction:
    - `estimate_salt_pairs(..., effective_volume_fraction=...)`
  - Summary now records `ion_effective_volume_fraction`.
- Updated workflow defaults:
  - `run_sim_1rkl.sh`: `BOX_PADDING_XY` default changed from `15.0` to `0.0` (keeps minimum-margin requirement without unnecessary extra XY growth).
- Validation:
  - Regenerated packed system:
    `python3 prepare_hybrid_system.py --output-dir outputs/hybrid_geom_fix_check2 ...`
  - Verified geometry from `outputs/hybrid_geom_fix_check2/hybrid_packed.MARTINI.pdb`:
    - protein spans: `x=34.612`, `y=11.844`
    - required lipid spans (>=3x protein span): `x>=103.836`, `y>=35.532`
    - actual lipid spans: `x=111.159`, `y=38.236` (passes)
    - lipid edges exactly match box edges in XY (passes).
  - Verified ion counts in summary:
    - `salt_pairs_target=13`, `NA=16`, `CL=13` (down from prior ~95/92 scale).
  - Syntax checks:
    - `python3 -m py_compile prepare_hybrid_system.py` passed.
    - `bash -n run_sim_1rkl.sh` passed.

## 2026-02-15 (Square XY Box + DOP Naming Canonicalization)
- Addressed new user constraints for generated `1rkl_hybrid.MARTINI.pdb`:
  - XY must be square (set shorter side to longer side).
  - Lipid naming must be `DOP` (not `DOPC`) in packed PDB output.
- Updated `prepare_hybrid_system.py`:
  - Added `canonical_lipid_resname(...)` and applied it during tiled lipid copy/output.
  - Changed target XY planning to a square footprint:
    - compute minimum required per-axis spans from protein extent + margin,
    - use `square_side = max(required_x, required_y)`,
    - crop target window to square centered on protein XY center.
  - Updated box assignment to enforce exact square XY at write time:
    - `set_box_from_lipid_xy(..., force_square_xy=True)` sets `box_x = box_y = max(lipid_span_x, lipid_span_y)`.
- Validation run:
  - `python3 prepare_hybrid_system.py ... --output-dir outputs/hybrid_square_check2`
  - Resulting PDB metrics:
    - `CRYST1`: `111.159 x 111.159 x 82.809` (XY square).
    - Residues: `DOP=3892`, `DOPC=0`.
    - Ions: `NA=41`, `CL=38` (`salt_pairs_target=38`, protein charge `-3`).
- Updated generated runtime file:
  - Copied `outputs/hybrid_square_check2/hybrid_packed.MARTINI.pdb` to `pdb/1rkl_hybrid.MARTINI.pdb`.

## 2026-02-15 (Minimum Box Z = 3x Protein Z Span)
- Added packed-box Z constraint in `prepare_hybrid_system.py`:
  - `set_box_from_lipid_xy(...)` now accepts `min_box_z`.
  - Main packing path computes `min_box_z_target = 3.0 * protein_z_span`.
  - Final `box_z` is set to `max(current_box_z, min_box_z_target)`.
- Added summary diagnostics:
  - `protein_z_span`
  - `min_box_z_target`
- Validation:
  - Generated `outputs/hybrid_square_z3_check/hybrid_packed.MARTINI.pdb`.
  - Verified `CRYST1`: `111.159 x 111.159 x 110.196`.
  - Verified protein Z span: `36.732`, so `3x = 110.196`.
  - Confirmed `box_z >= 3x protein_z_span` (exact equality in this run).
- Updated runtime file:
  - Copied `outputs/hybrid_square_z3_check/hybrid_packed.MARTINI.pdb` to `pdb/1rkl_hybrid.MARTINI.pdb`.

## 2026-02-15 (DOPC Name Preservation + 4-Character PDB Residue Handling)
- Addressed lipid-name mismatch where `DOPC` was being truncated to `DOP`:
  - `parse_pdb()` now reads residue names from `line[17:21]` (4-character support).
  - `write_pdb()` now writes residue names with 4-character width.
  - `canonical_lipid_resname()` now maps lipid inputs to `DOPC` (not `DOP`).
  - Lipid chain fallback remains blank in output (no forced `A`/`B` chain), preserving template-style `DOPC    <resid>` records.
- Regenerated packed output and replaced runtime file:
  - Source: `outputs/hybrid_square_z3_dopc_fmt_check/hybrid_packed.MARTINI.pdb`
  - Runtime target: `pdb/1rkl_hybrid.MARTINI.pdb`
- Verified key line format:
  - `ATOM     66  NC3 DOPC    1      ...`
- Verified geometry constraints still hold after formatting fix:
  - `CRYST1`: `111.159 x 111.159 x 110.196`
  - `box_z >= 3x protein_z_span` (exact in this run).

## 2026-02-15 (prepare_martini Residue Parse Compatibility Fix)
- Resolved stage preparation failure:
  - Error: `Unknown residue type 'AS' for atom 'SC1'` in `prepare_martini.py`.
  - Root cause: `prepare_martini.py` parsed residue name using 3-char slice (`line[17:20]`), which truncated protein residues when reading 4-char residue-field PDB formatting (`" ASN"` -> `"AS"`).
- Fix:
  - Updated residue parsing to use 4-character residue field:
    - `prepare_martini.py`: `residue_name = line[17:21].strip().upper()`
- Verification:
  - `python3 prepare_martini.py 1rkl_hybrid --stage minimization outputs/tmp_parse_fix` now succeeds.
  - Molecule summary confirms correct residue classification:
    - `PROTEIN: 1 chain(s) (31 residues)`
    - `DOPC: 278 molecules`
  - No unknown residue-type error after patch.

## 2026-02-15 (Bilayer Z-Centering + Fixed-in-Space Preproduction Hold)
- Updated packed-system Z placement in `prepare_hybrid_system.py`:
  - `set_box_from_lipid_xy(...)` now centers lipid bilayer midplane at `box_z/2`.
  - While centering, `box_z` is enlarged if needed so all atoms remain inside `[0, box_z]`.
  - Existing minimum `box_z >= 3 * protein_z_span` rule remains enforced.
  - Added summary diagnostics: `lipid_mid_z` and `box_half_z`.
- Regenerated runtime packed PDB:
  - `pdb/1rkl_hybrid.MARTINI.pdb` from `outputs/hybrid_zcenter_fixspace_check2/hybrid_packed.MARTINI.pdb`.
  - Verified: `lipid_mid_z == box_half_z == 55.098` (exact), `CRYST1 Z=110.196`.
- Updated C++ rigid-hold behavior for pre-production in `src/martini.cpp`:
  - Added per-engine reference coordinates for fixed atoms at fixation time.
  - `apply_fix_rigid_md(...)` now restores fixed atom positions to reference coordinates each step (in addition to zeroing forces/momenta).
  - This converts pre-production rigid hold from “free rigid drift” to “fixed in space”.
- Updated integration loop in `src/main.cpp`:
  - Re-applies `apply_fix_rigid_md(...)` after barostat/box updates each MD step so fixed atoms remain anchored even under NPT scaling.
- Build/verification:
  - Rebuilt with `cmake --build ../../obj -j4` successfully.
  - Runtime parse/log confirms hybrid preproduction fixed set is active (`preprod_fixed=343` in short check run).

## 2026-02-15 (Packing Gap Reduction Around Protein)
- Investigated user-observed protein–bilayer gap in `pdb/1rkl_hybrid.MARTINI.pdb`.
- Quantified root cause:
  - With previous default `--protein-lipid-cutoff 5.0`, nearest protein→lipid bead distance was ~`5.26 Å`, producing a visibly wide annular cavity.
  - Lipid removal logic removes whole residues if any bead is within cutoff; high cutoff amplifies cavity size.
- Adjustments made:
  - `prepare_hybrid_system.py` default `--protein-lipid-cutoff` changed `5.0 -> 3.0`.
  - `run_sim_1rkl.sh` default `PROTEIN_LIPID_CUTOFF` changed `5.0 -> 3.0`.
  - Protein initial centering now uses lipid-only coordinates (not ions) for translation reference.
- Validation after regeneration (`outputs/hybrid_gapfix_check/hybrid_packed.MARTINI.pdb`, copied to `pdb/1rkl_hybrid.MARTINI.pdb`):
  - nearest protein→lipid bead distance: `~3.40 Å` (from `~5.26 Å`)
  - median nearest distance: `~6.01 Å` (from `~8.67 Å`)
  - bilayer Z-centering remains exact (`lipid_mid_z == box_z/2`).

## 2026-02-16 (Stage 6.2 NaN Explosion Debug + Rigid Hold Stabilization)
- Investigated user-reported Stage 6.2 immediate blow-up (`potential nan`, rapid XY box growth).
- Verified script-level unit fix status:
  - `run_sim_1rkl.sh` already converts step counts to Upside time units (`duration = nsteps * dt`, `frame_interval = frame_steps * dt`).
  - User log (`--duration 500 --frame-interval 50`) was from an older invocation path.
- Reproduced current failure on the same checkpoint with corrected time units:
  - `outputs/martini_test_1rkl_hybrid/checkpoints/1rkl.stage_6.2.up` still reached `nan` early when hybrid preproduction rigid mode was enabled.
- Isolated fault domain by A/B control:
  - Same `.up` with `/input/hybrid_control.enable=0` remained stable (no NaN).
  - Confirms instability is in hybrid preproduction rigid handling, not baseline MARTINI stage setup.
- First mitigation tested:
  - Removed fixed-atom per-step coordinate rewrites from `src/martini.cpp` (`apply_fix_rigid_md` kept only force/momentum zeroing).
  - This removed NaNs but did not keep fixed atoms stationary (`max_fixed_disp ~2.0 Å` over stage-equivalent run).
- Final fix implemented:
  - Added integration-level fixed-atom freeze path in `src/deriv_engine.cpp`:
    - Build per-cycle fixed mask from `martini_fix_rigid::get_fixed_atoms(...)`.
    - In all integration-cycle variants (standard and mass-aware), skip momentum/position updates for fixed atoms and force momentum to zero.
  - Kept `src/martini.cpp::apply_fix_rigid_md` as force/momentum zeroing only (no hard position reset).
  - This enforces true fixed-in-space behavior without destabilizing coordinate rewrites.
- Rebuild:
  - `cmake --build ../../obj -j4` succeeded.
- Validation runs:
  - Short replay (`--duration 0.5 --frame-interval 0.1`) on copied Stage 6.2 checkpoint: stable, no NaN.
  - Stage-equivalent replay (`--duration 5 --frame-interval 0.5`): completed, no NaN, no barostat blow-up.
  - Fixed-atom displacement check on output trajectory:
    - `n_fixed=354`, `max_fixed_disp=0.0`, `mean_fixed_disp=0.0`.
- Files modified in this debug/fix pass:
  - `../../src/deriv_engine.cpp`
  - `../../src/martini.cpp`
  - `task_plan.md`

## 2026-02-16 (PO4 Z-Only Constraint + Protein Full Fix in Preproduction)
- Implemented requested preproduction constraint split:
  - Protein atoms remain fully fixed in space (`x,y,z` locked).
  - Lipid headgroup atoms (`PO4` roles) are constrained only along `z` (`x,y` free).
- Updated rigid-constraint registry in `../../src/martini.cpp`:
  - Added separate dynamic z-fixed atom set (`g_dynamic_z_fixed_atoms` / `g_z_fixed_atoms`).
  - Added API functions:
    - `set_dynamic_z_fixed_atoms(...)`
    - `clear_dynamic_z_fixed_atoms(...)`
    - `get_z_fixed_atoms(...)`
  - Updated hybrid preproduction atom selection:
    - `preprod_fixed_atom_indices` now includes only protein atoms.
    - new `preprod_z_fixed_atom_indices` includes non-protein `PO4` atoms.
  - Updated stage gating and cleanup to apply/clear both full-fixed and z-fixed sets.
  - Updated hybrid parse log to report both counts:
    - `preprod_fixed=<protein_count> preprod_zfixed=<po4_count>`.
  - Updated constraint application:
    - minimization: full-fixed zeroes all force components; z-fixed zeroes only `z` force.
    - MD: full-fixed zeroes all force/momentum; z-fixed zeroes only `z` force/momentum.
- Updated integrator enforcement in `../../src/deriv_engine.cpp`:
  - Added z-fixed mask retrieval from `martini_fix_rigid::get_z_fixed_atoms(...)`.
  - In all integration-cycle variants (mass-aware and unit-mass):
    - full-fixed atoms skip all momentum/position updates.
    - z-fixed atoms keep `x,y` updates, force `p_z=0`, and apply zero `z` displacement.
- Build:
  - `cmake --build ../../obj -j4` succeeded.
- Validation run:
  - Replayed Stage 6.2 checkpoint with hybrid preproduction active:
    - command used duration/time units (`--duration 1.0 --frame-interval 0.1`).
    - no NaN / no explosion; run completed.
    - runtime log confirmed split sets: `preprod_fixed=65 preprod_zfixed=289`.
  - Trajectory displacement checks from output H5:
    - protein `max_disp = 0.0` (fully fixed in space).
    - PO4 `max_abs_dz = 0.0` and `mean_abs_dz = 0.0` (z locked).
    - PO4 `max_dxy = 3.978` and `mean_dxy = 0.779` (`x,y` remain mobile).

## 2026-02-16 (Stage-Wise VTF Export + Mode 1/2 Support)
- Implemented extractor mode support in `extract_martini_vtf.py`:
  - Added CLI flag `--mode {1,2}`.
  - `mode 1`: writes all MARTINI particles to VTF/PDB output.
  - `mode 2`: writes non-protein MARTINI particles plus protein all-atom backbone (`N,CA,C,O`) rebuilt per BB from `input/hybrid_bb_map/reference_atom_coords` and live BB translations.
- Refactored extraction logic for robustness:
  - Added explicit handling for UPSIDE input/output position tensor layouts.
  - Added resilient box-length inference from `output/box`, log fallback, MARTINI attrs, and CRYST1 fallback.
  - Added bond remapping for filtered outputs and generated backbone bonds for mode 2 (`N-CA`, `CA-C`, `C-O`, and inter-residue `C-N`).
  - Kept both `.vtf` and `.pdb` output support.
- Updated workflow driver `run_sim_1rkl.sh`:
  - Added `extract_stage_vtf()` helper.
  - Generates VTF after every stage (`6.0`, `6.1`, `6.2`, `6.3`, `6.4`, `6.5`, `6.6`, `7.0`).
  - Uses mode policy requested:
    - stages `6.x` -> `--mode 1`
    - stage `7.0` -> `--mode 2`
  - Uses runtime packed PDB id (`$RUNTIME_PDB_ID`) for metadata alignment.
- Verification:
  - `python3 -m py_compile extract_martini_vtf.py` passed.
  - `bash -n run_sim_1rkl.sh` passed.
  - Functional extractor checks on existing stage file:
    - mode 1 run succeeded (`4216` output particles).
    - mode 2 run succeeded (`4275` output particles = non-protein MARTINI + backbone atoms).

## 2026-02-16 (Step-Count Duration Control; dt-Independent Stage Length)
- Implemented explicit step-count runtime control in `../../src/main.cpp`:
  - Added new CLI option: `--duration-steps <int>`.
  - `--duration-steps` now sets `n_round` directly (integration rounds), independent of `dt`.
  - Preserved backward compatibility: if `--duration-steps` is not provided, legacy `--duration` (time-based) behavior is unchanged.
  - `--duration` argument changed from required to optional; engine now requires either `--duration` or `--duration-steps`.
  - Progress display now reports steps when `--duration-steps` is used.
- Updated `run_sim_1rkl.sh` MD stage invocation:
  - Uses `--duration-steps "$nsteps"` instead of time-converted `--duration`.
  - Keeps frame interval handling unchanged.
  - Stage log now prints `duration(steps)=...` and `frame_steps=...`.
- Build/test:
  - Rebuilt binaries successfully: `cmake --build ../../obj -j4`.
  - Verified runtime behavior with a direct check:
    - command: `... --duration-steps 20 --time-step 0.010 ...`
    - output progress showed `0 / 20 ...`, confirming step-count mode is active and dt-independent for step count.

## 2026-02-16 (Audit: MARTINI Nonbonded Interaction Source = dry_martini_v2.1.itp)
- Performed source-path audit for MARTINI nonbonded setup:
  - `prepare_martini.py` loads nonbonded pair parameters from `ff_dry/dry_martini_v2.1.itp` via `read_martini3_nonbond_params(...)`.
  - Pair construction hard-fails on missing type-pair parameters (`ValueError`), so no silent fallback table is used.
  - Runtime `../../src/martini.cpp` consumes only HDF5 `pairs` + `coefficients` for LJ/Coulomb pair evaluation.
- Performed concrete coefficient consistency check on generated stage file:
  - target: `outputs/martini_test_1rkl_hybrid/checkpoints/1rkl.stage_6.2.up`
  - result: `pairs=8881396`, `missing_pairs=0`, `eps_sigma_mismatch=0` against `ff_dry/dry_martini_v2.1.itp`.
  - unique `(epsilon,sigma)` in `.up` all belonged to the dry table (`observed_not_in_expected=0`).
- Clarified Coulomb path:
  - LJ uses table-derived `epsilon/sigma`.
  - Coulomb is computed from pair charges (`q_i*q_j`) with the MARTINI Coulomb constant / Ewald controls in `../../src/martini.cpp`.

## 2026-02-16 (Fix: SC->BB Force Transfer Uses Martinize Bonded Topology Mapping)
- Audited hybrid sidechain force flow:
  - Probabilistic SC-environment pair forces use MARTINI pair coefficients (`epsilon/sigma/q`) from `martini_potential/coefficients`, sourced from `ff_dry/dry_martini_v2.1.itp`.
  - Runtime SC->BB transfer was previously residue-BB direct and did not prefer bonded-topology projection metadata.
- Implemented topology-driven transfer mapping in `prepare_hybrid_system.py`:
  - Added `--protein-itp` input.
  - Added martinize ITP parser for `[atoms]`, `[bonds]`, and `[constraints]` (constraints treated as very stiff edges).
  - Added CG<->ITP index mapping with sequence-order-first alignment to handle residue index offsets between packed PDB and ITP.
  - Added bonded-graph BB assignment for each SC proxy via shortest compliance path (`sum(1/k)`), using martinize bond force constants.
  - `hybrid_sc_map/proj_target_indices` and `proj_weights` now explicitly encode this bonded-topology transfer mapping (default single BB target, weight 1.0).
- Updated workflow wiring:
  - `run_sim_1rkl.sh` now passes `--protein-itp "${PROTEIN_ITP_EFFECTIVE}"` to `prepare_hybrid_system.py`.
- Updated runtime transfer preference in `../../src/martini.cpp`:
  - `project_sc_gradient_if_active(...)` now applies explicit `hybrid_sc_map` projection targets/weights first.
  - Residue-based BB lookup is retained as fallback only when explicit projection targets are unavailable.
- Validation:
  - `python3 -m py_compile prepare_hybrid_system.py` passed.
  - `bash -n run_sim_1rkl.sh` passed.
  - `cmake --build ../../obj -j4` passed.
  - Bonded BB mapping sanity check on martinize outputs:
    - `protein_cg_atoms=65`, `sc_nodes_mapped_to_bb_by_bonded_topology=34` (all SC proxies mapped).
  - End-to-end prep smoke test with new `--protein-itp` completed:
    - wrote `/tmp/hybrid_bonded_sc_map_check/hybrid_packed.MARTINI.pdb`
    - wrote `/tmp/hybrid_bonded_sc_map_check/hybrid_mapping.h5`
    - SC rows with missing primary BB target: `0`.

## 2026-02-16 (Stage Gating Audit: Backbone/Hybrid Activation Only in Production)
- Audited pre-production vs production activation path in C++ runtime:
  - `register_hybrid_for_engine(...)` computes `active = (current_stage == activation_stage)`.
  - Pre-production rigid mode applies full protein fixed-mask (`x,y,z`) while hybrid is inactive.
  - MARTINI potential hybrid projections (`refresh_bb_positions_if_active`, SC probabilistic projection, BB/SC gradient projection) are all guarded by `st.active`.
- Found and fixed a stage-transition edge case in minimization mode:
  - `../../src/main.cpp` previously forced `switch_simulation_stage(..., \"production\")` after minimization.
  - Updated behavior to restore the original pre-minimization stage instead.
  - This prevents accidental production-only activation during pre-production minimization runs.
- Added explicit runtime logging fields for verification:
  - `current_stage` and `hybrid_active` are now printed in the hybrid parse line in `../../src/martini.cpp`.
- Validation checks:
  - Build passed: `cmake --build ../../obj -j4`.
  - Pre-production file (`1rkl.stage_6.2.up`) reports:
    - `current_stage=minimization`, `activation_stage=production`, `hybrid_active=0`.
  - Production file (`1rkl.stage_7.0.up`) reports:
    - `current_stage=production`, `activation_stage=production`, `hybrid_active=1`.
  - Pre-production minimization run (`1rkl.stage_6.0.up --minimize`) now switches:
    - `minimization -> minimization` (restored), not to `production`.
  - Coordinate check on pre-production trajectory (`1rkl.stage_6.2.up`) confirms rigid hold:
    - protein `max_disp=0.0`, `mean_disp=0.0`.

## 2026-02-16 (Production Protein MARTINI Interaction Rules Enforced)
- Requirement checked and implemented for production-stage protein MARTINI interactions:
  - no BB-BB (bonded + nonbonded),
  - no SC-SC (bonded + nonbonded),
  - BB-SC only when both atoms belong to the same residue.
- Runtime implementation updates in `../../src/martini.cpp`:
  - Added atom role/residue metadata capture in hybrid runtime state (`atom_role_class`, `atom_residue_id`) from `/input/atom_roles` and `/input/residue_ids`.
  - Added pair-rule helpers for protein intra-interactions and multibody filtering.
  - Updated nonbonded pair skip logic (`skip_pair_if_intra_protein`) to apply selective rule instead of global intra-protein exclusion.
  - Applied same rule to bonded/multibody terms during hybrid-active production:
    - `dist_spring` filters disallowed protein pairs.
    - `angle_spring` and `dihedral_spring` filter any term containing disallowed protein atom-pair combinations.
- Build/test:
  - Rebuilt successfully: `cmake --build ../../obj -j4`.
  - Data-level audit on production file `outputs/martini_test_1rkl_hybrid/checkpoints/1rkl.stage_7.0.up`:
    - protein bonded pairs in file:
      - BB-BB: `30`, SC-SC: `9`, BB-SC same-res: `28`.
      - runtime rule keeps only BB-SC same-res; filters BB-BB + SC-SC.
    - protein intra nonbonded pairs in file:
      - total `2013`, BB-BB `435`, SC-SC `552`, BB-SC same-res `6`, BB-SC diff-res `1020`.
      - runtime rule allows only `6` same-res BB-SC pairs; all others filtered.

## 2026-02-16 (Utility: Production-Only Stage 7 Runner for Faster Debug)
- Added `test_prod_run_sim_1rkl.sh` to run only production stage `7.0`, assuming stages `6.0-6.6` and hybrid-prep artifacts already exist.
- Script behavior mirrors the production block in `run_sim_1rkl.sh`:
  - validates required artifacts (`stage_6.6.up`, runtime hybrid PDB/ITP, `hybrid_mapping.h5`),
  - prepares `stage_7.0.prepared.up` with production stage label and hybrid mapping injection,
  - initializes `stage_7.0.up` from `stage_6.6.up`,
  - runs production MD using `--duration-steps` and time-based frame interval conversion,
  - extracts production VTF in mode `2`.
- File permissions and syntax checks:
  - `chmod +x test_prod_run_sim_1rkl.sh`
  - `bash -n test_prod_run_sim_1rkl.sh` passed.

## 2026-02-17 (Rotamer Numbering Consistency Audit)
- Re-read `task_plan.md` and audited rotamer-state encoding/decoding paths across Python and C++:
  - `../../py/upside_config.py` encodes `id_seq` as `(node_id << 8) | (n_rot << 4) | rot` with `n_bit_rotamer=4`.
  - `../../src/rotamer.cpp` decodes IDs with the same bit layout and exposes runtime logs `node_lookup=[n_rotamer,node_id]` and `node_marginal` (6-wide rows, active states in `[0, n_rot)`).
  - `../../src/martini.cpp` hybrid decode and probabilistic refresh use the same `(rot, n_rot, node_id)` convention.
- Ran runtime cross-check on `../01.GettingStarted/inputs/1UBQ.up`:
  - `placement_fixed_point_vector_only` rows: `296`
  - decoded invalid rows: `0`
  - unique placement keys `(n_rot,node_id)`: `76`
  - rotamer `n_node`: `76`
  - key mismatch counts: `placement_missing_in_lookup=0`, `lookup_missing_in_placement=0`
  - row-level validity: `ok=296`, `bad=0`
- Verified sidechain-library state indexing used by `predict_chi1.py`:
  - `../../parameters/ff_2.1/sidechain.h5` `restype_and_chi_and_state[:,3]` spans `0..5` and starts at `0` for every restype (no 1-based offset).
  - `predict_chi1.py` uses these state IDs directly as array indices (`probs[s]`), consistent with 0-based rotamer numbering.
- Confirmed current `1rkl` stage checkpoint status:
  - `outputs/martini_test_1rkl_hybrid/checkpoints/1rkl.stage_6.2.up` has no `rotamer` or `placement*_point_vector_only` nodes.
  - `hybrid_sc_map/rotamer_id` is single-state (`min=max=0`, all probabilities `1.0`), so live probabilistic rotamer weighting is not active in that file.

## 2026-02-17 (Strict No-Fallback SC Mapping + Rigid Rotamer Transform)
- Applied user-requested behavior in `../../src/martini.cpp`:
  - Removed fallback SC placement/probability behavior in production hybrid coupling.
  - Added strict production checks requiring both `rotamer` and `placement*_point_vector_only` runtime nodes for probabilistic SC mapping.
  - Ensured placement groups are rebuilt after runtime node discovery (previously they were parsed before node-name discovery, yielding empty groups).
  - Enabled SC row expansion for multi-bead residues (removed one-row-per-residue limitation in placement expansion gate).
- Implemented rigid per-rotamer SC placement:
  - Added residue reference placement-group selection (prefer rotamer `0`, else first group).
  - Added rigid transform mapping from reference placement frame to current rotamer placement frame.
  - SC row local offsets are initialized once from current proxy geometry when `hybrid_sc_map/local_pos` is zero, then transformed rigidly each step.
  - SC position build now uses placement-only mapping; static `hybrid_sc_map` target-position fallback is removed.
- Tightened probabilistic rotamer weighting:
  - `refresh_sc_row_probabilities_from_rotamer(...)` now returns success/failure.
  - MARTINI compute path now aborts if live rotamer probabilities cannot be resolved for every SC row.
- Updated `task_plan.md` to record the revised strict behavior and new blocker implication for old stage files lacking rotamer/placement nodes.
- Build and validation:
  - Rebuilt successfully: `cmake --build ../../obj -j4`.
  - Runtime sanity check on current production file:
    - Command: `../../obj/upside outputs/martini_test_1rkl_hybrid/checkpoints/1rkl.stage_7.0.up --duration-steps 0 ...`
    - Result: expected fast-fail with `ERROR: Hybrid production SC coupling requires placement*_point_vector_only node`.

## 2026-02-17 (Production Rotamer-Node Injection Bug Fix in Run Scripts)
- Fixed production-stage augmentation logic in both workflow scripts:
  - `run_sim_1rkl.sh`
  - `test_prod_run_sim_1rkl.sh`
- Root cause:
  - augmentation code assumed `sidechain.h5` contains `rotamer_prob_fixed`, but current `ff_2.1/sidechain.h5` provides `rotamer_prob` only.
  - this raised `KeyError` during stage-7 preparation, so production files were not augmented with rotamer/placement nodes.
- Implemented fix:
  - if `rotamer_prob_fixed` exists, use it directly as fixed scalar placement energy.
  - otherwise derive fixed scalar energy from `rotamer_prob` by phi/psi averaging and `-log(mean_prob)` with clipping for numerical safety.
- Validation:
  - `bash -n run_sim_1rkl.sh` passed.
  - `bash -n test_prod_run_sim_1rkl.sh` passed.
  - `PROD_70_NSTEPS=1 PROD_FRAME_STEPS=1 ./test_prod_run_sim_1rkl.sh` now injects nodes successfully:
    - log shows `Injected production rotamer nodes ...`
    - runtime log shows `placement_node=placement_fixed_point_vector_only rotamer_node=rotamer`.
    - original error `Hybrid production SC coupling requires placement*_point_vector_only node` no longer occurs.
  - Additional direct init check passed:
    - `../../obj/upside ...stage_7.0.up --duration-steps 0 ...` exits cleanly with hybrid production parsing active.
- Remaining issue observed (separate from missing-node bug):
  - 1-step production run still crashes with `Bus error: 10` after entering integration; this indicates a runtime stability/memory issue after initialization, not a stage-7 node-injection failure.

## 2026-02-17 (Bus Error Debug + Fix for Stage 7 Production)
- Investigated crash report `~/Library/Logs/DiagnosticReports/upside-2026-02-17-003758.ips`.
- Crash stack showed SIGBUS in rotamer spline evaluation path:
  - `bead_interaction.h` -> `spline.h::deBoor_value_and_deriv` -> `rotamer.cpp::fill_holders`.
- Root cause:
  - NaN/invalid SC placement coordinates were generated from degenerate `affine_alignment` triplets in script-generated production nodes.
  - Two residues had duplicate triplet endpoints (`a == c`), producing zero-area reference frames and invalid rotamer geometry.
- Implemented fix in both scripts:
  - `run_sim_1rkl.sh`
  - `test_prod_run_sim_1rkl.sh`
  - Replaced simple neighbor triplet selection with robust triplet construction:
    - build candidate atoms from nearby BBs + SC proxy + global fallbacks,
    - enforce distinct atoms (`a != bb != c`, `a != c`),
    - choose the pair maximizing triangle area around BB to avoid near-collinear/degenerate frames.
- Validation:
  - `bash -n run_sim_1rkl.sh` passed.
  - `bash -n test_prod_run_sim_1rkl.sh` passed.
  - `PROD_70_NSTEPS=1 PROD_FRAME_STEPS=1 ./test_prod_run_sim_1rkl.sh` completed without SIGBUS.
  - `PROD_70_NSTEPS=10 PROD_FRAME_STEPS=10 ./test_prod_run_sim_1rkl.sh` completed without SIGBUS.
- Current status:
  - Bus error is resolved.
  - Production still shows rapid energy blow-up / NaN potential (physical stability issue), but no immediate runtime memory crash.

## 2026-02-17 (Rg Status-Line Fix: Use All-Atom Backbone Mapping)
- Updated runtime status-line `Rg` calculation in `../../src/main.cpp` to avoid legacy `3*n_res` first-atom assumption.
- Implemented per-system backbone atom index collection at initialization:
  - primary source: `/input/hybrid_bb_map/atom_indices` + `/input/hybrid_bb_map/atom_mask` (captures mapped `N,CA,C,O` all-atom backbone indices),
  - fallback: legacy `3*n_res` sequence-based indexing when hybrid mapping is absent.
- Added cached `System::rg_backbone_atom_indices` and reused it during frame printout for `Rg` computation.
- Validation:
  - Rebuilt successfully: `cmake --build ../../obj -j4`.
  - Production smoke run (`PROD_70_NSTEPS=1`) now prints numeric Rg:
    - `... hbonds, Rg  12.6 A, potential ...`
  - `Rg` no longer prints `N/A` in this hybrid workflow.

## 2026-02-17 (SC-Environment Stabilization: Force Caps + 200-Step Inner Relaxation)
- Implemented requested production hybrid SC-environment stabilization in `../../src/martini.cpp`:
  - Added SC-env force caps (LJ and Coulomb) applied specifically in probabilistic SC-environment interactions.
  - Replaced single-evaluation SC-environment force path with an inner SC-only relaxation loop:
    - environment and backbone coordinates are fixed during inner loop,
    - SC rows are restrained near mapped rotamer positions,
    - inner relaxation runs `200` steps by default,
    - only final relaxed-step weighted SC force is projected back to BB/environment.
- Added hybrid-control runtime parameters (with defaults) to `HybridRuntimeState` and parser:
  - `sc_env_lj_force_cap` (default `25.0`)
  - `sc_env_coul_force_cap` (default `25.0`)
  - `sc_env_relax_steps` (default `200`)
  - `sc_env_relax_dt` (default `0.002`)
  - `sc_env_restraint_k` (default `5.0`)
  - `sc_env_max_displacement` (default `2.0`)
- Runtime startup log now prints these SC controls for auditability.
- Validation:
  - Rebuilt successfully: `cmake --build ../../obj -j4`.
  - `PROD_70_NSTEPS=2` production smoke run:
    - shows new controls in log (`sc_cap_lj`, `sc_cap_coul`, `sc_relax_steps=200`, etc.),
    - no crash/NaN in the tested window.
  - `PROD_70_NSTEPS=10` check:
    - no immediate NaN explosion observed over 10 steps,
    - MARTINI potential remains very large (physical calibration still needed).

## 2026-02-17 (Audit: Pre-Production Rigidity + Potential Output Semantics)
- Traced workflow and runtime paths relevant to the user report:
  - `run_sim_1rkl.sh` stage labels and extraction modes (`6.x -> minimization`, `7.0 -> production`, VTF mode `1` vs `2`).
  - `../../src/main.cpp` logger wiring for `output/potential`.
  - `../../src/deriv_engine.cpp` potential accumulation behavior.
  - `../../src/martini.cpp` hybrid pre-production rigid-mask construction and stage gating.
- Runtime artifact checks on `outputs/martini_test_1rkl_hybrid/checkpoints` and logs:
  - Stage logs 6.0-6.6 report:
    - `current_stage=minimization`
    - `activation_stage=production`
    - `hybrid_active=0`
    - `preprod_mode=rigid`
    - `preprod_fixed=65` (protein) and `preprod_zfixed=289` (`PO4`)
  - Stage files 6.0-6.6 contain scalar `output/potential` datasets with shape `(n_frame, 1)` and no protein-only potential channel.
- Quantitative rigidity verification (`h5py` checks):
  - Protein displacement is exactly zero in all pre-production stages:
    - stage 6.0-6.6: `max |r - r_ref| = 0.000000e+00 Å`
    - cross-stage continuity 6.0->6.6 for protein: `max |stage_start - prev_stage_end| = 0.000000e+00 Å`
  - `PO4` headgroup behavior matches z-only constraint:
    - `max |dz| = 0.000000e+00 Å` while `|dx|, |dy|` are non-zero.
- Conclusion recorded:
  - Pre-production rigid hold is active and effective throughout stages 6.0-6.6.
  - Reported potential drift is expected because `output/potential` is total system potential, not protein all-atom-backbone Upside-only potential.

## 2026-02-17 (Potential Split + Pre-Production NPT/Pressure Fixes)
- Implemented split potential logging in `../../src/main.cpp`:
  - Added runtime split of potential terms into MARTINI vs non-MARTINI buckets by node prefix.
  - Changed `output/potential` logger to record non-MARTINI (Upside-side) contribution.
  - Added new `output/martini_potential` logger for MARTINI contribution.
  - Updated verbose console print to include all three values:
    - `potential` (Upside-side),
    - `martini_potential`,
    - `total` (engine sum).
- Implemented pre-production NPT compatibility with rigid hold in `../../src/box.cpp`:
  - Removed hard skip that previously disabled barostat updates when pre-production rigid mode was active.
  - Updated semi-isotropic scaling to honor fixed masks:
    - fully fixed atoms are not scaled,
    - z-fixed atoms scale in `x/y` only (z unchanged).
- Updated stage NPT targets in `run_sim_1rkl.sh` to follow CHARMM-GUI stage intent with explicit unit-converted values:
  - Added `BAR_1_TO_EUP_PER_A3=0.000020659477`.
  - Added `set_stage_npt_targets()` and wired calls before each stage 6.0-6.6 preparation:
    - 6.0-6.1 use `0 bar`,
    - 6.2-6.6 use `1 bar` (`semi-isotropic`, `compressibility_z=0` retained).
- Validation run:
  - Rebuilt binaries: `cmake --build ../../obj -j4` succeeded.
  - Probe run on copied stage file:
    - `../../obj/upside /tmp/1rkl.stage_6.2.npt_probe.up --duration-steps 100 ...`
    - runtime output confirms split print and NPT update line:
      - `Initial potential energy (Upside/MARTINI/Total): 0.00/290333.76/290333.75`
      - `[NPT] t 0.500 box ...`
  - H5 output verification on probe file:
    - `output/martini_potential` dataset present.
    - `output/potential` present and zero for this MARTINI-only stage (expected).
    - `output/pressure` non-zero at later frame (`[0.00385833, -0.00927570]`).
    - protein remained rigid (`max displacement = 0.0 A`).
    - PO4 remained z-fixed (`max |dz| = 0.0 A`, non-zero `|dx|/|dy|`).
  - Additional production-initialization probe (`--duration-steps 0`) on stage-7 copy:
    - `Initial potential energy (Upside/MARTINI/Total): 177.63/153485216298.49/153485213696.00`
    - confirms non-MARTINI bucket is populated when production rotamer/placement nodes are active.

## 2026-02-18 (Preparation Patch: Explicit SC Control Attrs)
- Patched production-stage preparation in both workflow scripts:
  - `run_sim_1rkl.sh`
  - `test_prod_run_sim_1rkl.sh`
- Added script-level configuration variables for SC coupling control attrs:
  - `SC_ENV_LJ_FORCE_CAP` (default `25.0`)
  - `SC_ENV_COUL_FORCE_CAP` (default `25.0`)
  - `SC_ENV_RELAX_STEPS` (default `200`)
  - `SC_ENV_RELAX_DT` (default `0.002`)
  - `SC_ENV_RESTRAINT_K` (default `5.0`)
  - `SC_ENV_MAX_DISPLACEMENT` (default `2.0`)
- Added `set_hybrid_sc_controls()` helper to both scripts to write these values into `/input/hybrid_control` attrs in the prepared stage file.
- Wired helper into `prepare_stage_file()` for production stage so attr injection happens during preparation, before production MD runs.
- Validation:
  - `bash -n run_sim_1rkl.sh` passed.
  - `bash -n test_prod_run_sim_1rkl.sh` passed.
  - `rg` checks confirm function presence and production-stage call sites in both scripts.
  - Runtime preparation probe:
    - `PROD_70_NSTEPS=0 PROD_FRAME_STEPS=1 bash test_prod_run_sim_1rkl.sh` completed.
    - Stage-7 startup log reports explicit SC controls in parsed hybrid config:
      - `sc_cap_lj=25.000`
      - `sc_cap_coul=25.000`
      - `sc_relax_steps=200`
      - `sc_relax_dt=0.0020`
      - `sc_rest_k=5.000`
      - `sc_max_disp=2.000`
  - HDF5 attribute check on regenerated `outputs/martini_test_1rkl_hybrid/checkpoints/1rkl.stage_7.0.up` confirms all six SC control attrs are present under `/input/hybrid_control`.

## 2026-02-19 (Production Instability Root-Cause Audit)
- Investigated reported production-stage behavior in `outputs/martini_test_1rkl_hybrid/checkpoints/1rkl.stage_7.0.up`:
  - large MARTINI potential spikes,
  - rapid protein-backbone disruption (Rg growth).
- Artifact/structure checks:
  - `n_atom=4216`, protein membership count `65` (MARTINI protein beads only), environment `4151`.
  - Potential node set is MARTINI-centric (`dist_spring`, `angle_spring`, `dihedral_spring`, `martini_potential`) plus injected rotamer/placement nodes; no all-atom Upside backbone node family is present in this stage artifact.
  - `hybrid_bb_map` and `hybrid_sc_map` projection entries are MARTINI-index mapped with weights `[1,0,0,0]` (no direct all-atom `N/CA/C/O` distribution in this artifact).
- Trajectory diagnostics from stage-7 output frames:
  - Protein Rg grows monotonically (`~12.95 Å` at step 0 to `~28.50 Å` at step 4500).
  - Protein bond geometry strongly diverges from equilibrium:
    - mean `|d - d0|` rises from `0.223 Å` (step 0) to `19.401 Å` (step 4500),
    - max `|d - d0|` reaches `73.838 Å`.
- C++ code-path root causes identified:
  - Production hybrid filtering (`allow_protein_pair_by_rule`) disallows `BB-BB`; this rule is used not only in nonbonded MARTINI but also in bonded/multibody potentials:
    - `DistSpring` skips pair if `allow_intra_protein_pair_if_active` fails,
    - `AngleSpring`/`DihedralSpring` skip if `allow_multibody_term_if_active` fails.
  - Net effect in current artifact: protein bonded constraints are skipped in production while no all-atom backbone force field is actively replacing them.
  - SC force capping is force-only, not potential-only:
    - `eval_pair_force` caps force vectors but still accumulates raw LJ/Coulomb potential, so short-distance SC-environment contacts can create large instantaneous potential spikes.
- Quantitative spike support (pair scan on saved frames):
  - protein-environment nearest contacts reach `~1.29–2.07 Å`,
  - corresponding uncapped LJ pair potentials can be very large (`~1.7e7` at step 1500, `~1.2e5` at step 2000, before probabilistic weighting).
  - This is consistent with observed total MARTINI potential volatility when weighted SC terms vary by rotamer probabilities.

## 2026-02-19 (BB Mapping Audit: Reference AA vs Active Projection)
- Re-checked preparation and runtime BB mapping flow against requested behavior.
- Preparation (`prepare_hybrid_system.py`) currently writes:
  - reference all-atom backbone metadata under `hybrid_bb_map`:
    - `reference_atom_indices` (`N/CA/C/O`)
    - `reference_atom_coords`
    - `bb_comment`
  - active projection mapping for runtime as MARTINI self-mapping:
    - `atom_indices=[bb_atom_index,-1,-1,-1]`
    - `weights=[1,0,0,0]`
- Runtime (`src/martini.cpp`) consumes only active projection datasets (`atom_indices`, `atom_mask`, `weights`, `bb_atom_index`) for COM refresh and BB-force projection; reference AA metadata is not used in force transfer.
- Verified in generated production artifact `outputs/martini_test_1rkl_hybrid/checkpoints/1rkl.stage_7.0.up`:
  - `reference_atom_indices` and `reference_atom_coords` are present.
  - active `atom_indices` rows are MARTINI BB self-indices, and `weights` are `[1,0,0,0]`.
- Net current behavior:
  - BB force is projected back onto BB itself (self-target), not distributed to all-atom `N/CA/C/O` components.
  - Requested “block of 4 components” force splitting is therefore not active in this artifact path.

## 2026-02-19 (BB Force-Split Activation, No Direct-BB Projection)
- Implemented preparation/runtime changes to disable direct BB force application and activate component-wise BB force transfer.
- `prepare_hybrid_system.py` updates:
  - Added BB component definitions (`N/CA/C/O` + mass ratios).
  - `collect_bb_map()` now emits active `hybrid_bb_map/atom_indices` + `weights` for four component targets with normalized mass-fraction weights.
  - Direct BB self-target entries (`atom_indices == bb_atom_index`) are excluded.
  - Added nearest non-BB protein fallback selection to keep four active component targets when direct all-atom indices are not representable in current runtime index space.
  - `bb_comment` now records both reference all-atom indices and active runtime targets/weights.
- `src/martini.cpp` updates:
  - `refresh_bb_positions_if_active()` skips BB self-target terms when rebuilding BB from component targets.
  - `project_bb_gradient_if_active()` skips BB self-target terms when distributing BB gradient to component targets.
  - Net: no BB self-force transfer path remains, even if legacy mappings still contain BB self-target terms.
- Verification:
  - Python compile check passed: `python3 -m py_compile prepare_hybrid_system.py`.
  - C++ rebuild passed: `cmake --build ../../obj -j4`.
  - Regenerated mapping validates: `validate_hybrid_mapping.py ... --n-atom 4216` -> OK.
  - Regenerated mapping statistics:
    - `n_bb=31`
    - `rows_with_any_active_target=31`
    - `rows_with_4_active_targets=31`
    - `bb_self_target_terms=0`
  - Activated in workflow output by regenerating `outputs/martini_test_1rkl_hybrid/hybrid_prep/hybrid_mapping.h5` and re-running `test_prod_run_sim_1rkl.sh` (`PROD_70_NSTEPS=0`), then confirming stage file `outputs/martini_test_1rkl_hybrid/checkpoints/1rkl.stage_7.0.up` has:
    - 31 BB rows, all 4-way active mappings,
    - zero BB self-target terms.

## 2026-02-18 (Production Instability Isolation: 1-Step Stage-7 Run)
- Reproduced requested one-step production run:
  - `PROD_70_NSTEPS=1 PROD_FRAME_STEPS=1 ./test_prod_run_sim_1rkl.sh`
  - Observed startup energies:
    - `Initial potential energy (Upside/MARTINI/Total): 177.63/2418077075941.53/2418077073408.00`
  - Confirms instability is already present at step `0` and dominated by `martini_potential`.
- Performed controlled toggles on copied stage-7 files (`--duration-steps 0`) to isolate source:
  - Hybrid disabled (`/input/hybrid_control/enable=0`) -> MARTINI potential `-14740.32`.
  - Hybrid enabled but `hybrid_sc_map/rotamer_id` set empty (`n_sc=0`) -> MARTINI potential `3155.49`.
  - Hybrid enabled with SC mapping active (`n_sc=153`) -> MARTINI potential `2418077075941.53`.
  - Conclusion: blow-up originates in production SC probabilistic coupling path, not generic MARTINI nonbonded or hard-sphere setting.
- Ran SC-row ablation sweep (single-row `hybrid_sc_map` runs, `--duration-steps 0`) to identify dominant offender:
  - Worst row: `row=6` (`residue_index=10`, `proxy_atom_index=13`, LEU SC1), `n_sc=3`, MARTINI potential `2417895410149.53`.
  - Next largest rows are far smaller (`1.53e8`, `1.43e7`, ...), indicating one dominant pathological row.
- Verified single-row dominance:
  - Forcing only row `6` rotamer `2` (`n_sc=1`) still yields MARTINI potential `2417945479653.53`.
- Root-cause analysis in C++ mapping math (`../../src/martini.cpp`):
  - `build_sc_row_proxy_pos` uses `local = proxy - ref_centroid`, computes `t = tgt_centroid - R*ref_centroid`, then applies `mapped = R*local + t`.
  - This is inconsistent for centroid-relative `local` coordinates (target centroid is not restored), creating large translation errors in mapped SC positions.
  - For the dominant LEU row, this produces a near-overlap with lipid atom pair `(13,2050)` (`eps=0.92625856`, `sig=4.7`), with reconstructed distance ~`0.1715 Å` to DOPC `D3B`, which is sufficient to produce ~`1e12`-scale LJ repulsion.

## 2026-02-18 (Fix Applied: SC Centroid Restore in C++)
- Patched `../../src/martini.cpp` in `build_sc_row_proxy_pos(...)`:
  - changed SC mapping from `mapped = R*local + t` to `mapped = R*local + tgt_centroid`.
  - this matches `local` semantics (`proxy - ref_centroid`) and correctly restores target-frame translation.
- Rebuilt binaries after patch:
  - `source ../../.venv/bin/activate && source ../../source.sh && cmake --build ../../obj -j4`
  - Build succeeded for `upside`, `upside_calculation`, `upside_engine` (warnings unchanged).
- Re-ran 1-step production:
  - `PROD_70_NSTEPS=1 PROD_FRAME_STEPS=1 ./test_prod_run_sim_1rkl.sh`
  - Initial energies changed from:
    - before fix: MARTINI `2418077075941.53`
    - after fix: MARTINI `451819211237.53`
  - indicates substantial reduction (~5.35x lower), but instability remains.
- Post-fix targeted checks:
  - `n_sc=0` baseline remains stable: MARTINI `3155.49`.
  - Former dominant row (`row=6`, `residue_index=10`, `proxy_atom_index=13`) is resolved:
    - before fix (single-row): `2417895410149.53`
    - after fix (single-row): `3155.02`.
  - New dominant row after fix:
    - `row=26` (`residue_index=26`, `proxy_atom_index=49`) single-row MARTINI `451816163813.53`.

## 2026-02-18 (Follow-up Fix: SC Row-Expansion Capacity Guard)
- Continued debug for remaining dominant row (`residue 26 / proxy 49`) showed a representation mismatch:
  - placement groups for residue `26` carry only `1` placement point per rotamer state,
  - but `hybrid_sc_map` base rows include `3` MARTINI proxies for that residue.
- Control experiment on copied stage file:
  - dropping residue-26 extra proxies (`48`, `49`) while keeping proxy `47` reduced step-0 MARTINI potential to `221225.53`.
- Implemented runtime guard in `../../src/martini.cpp`:
  - added `compute_sc_proxy_limit_from_placement(...)`.
  - in `expand_sc_rows_from_placement(...)`, per-residue base proxy rows are capped by placement-point capacity per rotamer (minimum group size across placement groups).
  - excess base proxy rows are skipped before rotamer-state expansion.
- Rebuilt after code changes:
  - `source ../../.venv/bin/activate && source ../../source.sh && cmake --build ../../obj -j4`
  - build succeeded (`upside`, `upside_calculation`, `upside_engine`; existing warnings unchanged).
- Validation after follow-up fix:
  - `PROD_70_NSTEPS=1 PROD_FRAME_STEPS=1 ./test_prod_run_sim_1rkl.sh`
    - runtime reports `n_sc=117` (reduced from `153`).
    - step-0 MARTINI potential: `219631.41` (no catastrophic startup blow-up).
  - `PROD_70_NSTEPS=5 PROD_FRAME_STEPS=1 ./test_prod_run_sim_1rkl.sh`
    - run completed, all energies finite.
    - MARTINI potential range observed in status lines: `1.16e5` to `5.46e6`.
  - direct continuation probe:
    - `../../obj/upside ... --duration-steps 20 ...`
    - run completed, no NaN/crash.
    - MARTINI potential remained finite; observed spikes up to `2.45e7` (still elevated but far below prior `1e11-1e12` startup regime).

## 2026-02-18 (Continued Debug: Residual Spike Isolation + Additional Runtime Fixes)
- Performed frame-accurate spike analysis using frozen-step-0 `local_pos` replay:
  - identified that naive frame replay (`--duration-steps 0`) re-initializes `local_pos` and can hide true in-run spike states.
  - reproduced spike frame potentials by forcing `local_pos` to step-0-derived values.
- Leave-one-out ablation on spike frames isolated dominant contributor:
  - frame `13` (`martini ~2.37e7`): dropping `row 15` (`residue_index=19`, `proxy_atom_index=31`) reduced MARTINI to `~5.16e3`.
  - frame `17` (`martini ~1.44e9` in later variant): dropping the same row reduced MARTINI to `~5.14e4`.
- Implemented additional runtime fixes in `../../src/martini.cpp`:
  - Reworked SC proxy-capacity handling:
    - preserved all expanded SC rows,
    - added per-residue capacity map (`sc_proxy_limit_by_residue`),
    - added per-step proxy subset selection by nearest-environment-distance safety metric when proxies exceed representable capacity.
  - Added underdetermined-placement handling (`<3` reference points):
    - refresh `sc_local_pos` every step from current proxy coordinates (instead of relying on stale one-time initialization),
    - cap mapped displacement from current proxy by `sc_env_max_displacement` during SC target construction.
- Build and validation:
  - rebuilt after each patch: `cmake --build ../../obj -j4` (success).
  - intermediate 20-step run after first added fix removed earlier mid-run `~2.4e7` spikes but produced a late spike (`~1.44e9`).
  - after underdetermined local-refresh + displacement-cap patch:
    - `PROD_70_NSTEPS=20 PROD_FRAME_STEPS=1 ./test_prod_run_sim_1rkl.sh` completed with no NaN/crash,
    - MARTINI potential stayed in a low finite band (`~8e3` to `~5.3e4`), with prior `1e7-1e9` spike regimes absent in this probe.

## 2026-02-19 (PDB N/CA/C/O Index Preservation + Runtime Carrier Mapping)
- Re-checked user concern: stage artifacts used only MARTINI protein beads (`protein_membership=65`), so direct PDB `N/CA/C/O` indices were not usable as active runtime targets.
- Updated `prepare_hybrid_system.py`:
  - `collect_bb_map()` now stores active BB component indices in `protein_aa_pdb_0based` index space (no nearest-MARTINI fallback).
  - Added `hybrid_bb_map` attrs:
    - `atom_index_space=protein_aa_pdb_0based`
    - `runtime_index_space=stage_runtime_after_injection`
  - Applied the same translation shift to `protein_aa_atoms` as to MARTINI protein coordinates so `reference_atom_coords` are in the packed frame.
- Updated `validate_hybrid_mapping.py`:
  - Added BB index-space decoding.
  - In `protein_aa_pdb_0based` mode, BB target validation uses non-negative reference-index checks (runtime bounds checks are deferred to stage-file injection/runtime mapping).
  - Validator output now reports `bb_index_space`.
- Reworked `inject_hybrid_mapping()` in both workflow scripts:
  - `run_sim_1rkl.sh`
  - `test_prod_run_sim_1rkl.sh`
  - New behavior:
    - copy hybrid groups as before;
    - append AA reference-coordinate carrier particles to stage `/input` arrays using preserved PDB index space (`runtime_index = base_n_atom + pdb_index`);
    - extend required per-atom datasets (`pos`, `mom`, `vel`, `mass`, `charges`, `type`, `atom_names`, `atom_roles`, `residue_ids`, `molecule_ids`, `particle_class`);
    - extend `/input/potential/martini_potential/{atom_indices,charges}` lengths to match augmented `n_atom`;
    - rewrite `/input/hybrid_bb_map/{atom_indices,atom_mask,weights}` to runtime carrier indices and mass-fraction weights;
    - rewrite `/input/hybrid_env_topology/protein_membership` to augmented length and mark carriers as protein.
- Patched `set_initial_position.py` to handle source/target atom-count mismatch:
  - on mismatch, copy overlap and preserve target extra coordinates (prevents losing injected AA carriers when seeding from older checkpoints).
- Validation and smoke checks:
  - `python3 -m py_compile prepare_hybrid_system.py validate_hybrid_mapping.py set_initial_position.py` passed.
  - `bash -n run_sim_1rkl.sh` and `bash -n test_prod_run_sim_1rkl.sh` passed.
  - Generated and validated mapping:
    - `python3 prepare_hybrid_system.py ... --output-dir outputs/hybrid_indexspace_check`
    - `python3 validate_hybrid_mapping.py outputs/hybrid_indexspace_check/hybrid_mapping.h5`
    - result reports `bb_index_space=protein_aa_pdb_0based`.
  - Production-only prep smoke:
    - `PROD_70_NSTEPS=0 PROD_FRAME_STEPS=1 ./test_prod_run_sim_1rkl.sh` now succeeds.
    - Stage file confirms:
      - `n_atom=4449` (augmented from 4216),
      - `reference_index_offset=4216`, `reference_index_count=233`,
      - zero BB self-target terms,
      - BB active mapping satisfies `atom_indices = offset + reference_atom_indices`.
  - Reduced full-workflow smoke:
    - `MIN_60_MAX_ITER=1 MIN_61_MAX_ITER=1 EQ_62_NSTEPS=0 EQ_63_NSTEPS=0 EQ_64_NSTEPS=0 EQ_65_NSTEPS=0 EQ_66_NSTEPS=0 PROD_70_NSTEPS=0 EQ_FRAME_STEPS=1 PROD_FRAME_STEPS=1 ./run_sim_1rkl.sh`
    - completed through stage `7.0` with augmented `n_atom=4449` across checkpoints.
  - Functional production-step check:
    - `PROD_70_NSTEPS=1 PROD_FRAME_STEPS=1 ./test_prod_run_sim_1rkl.sh`
    - production stage executed for one step with hybrid active and augmented atom count (`n_atom=4449`) without initialization/runtime schema errors.

## 2026-02-19 (Production Upside Backbone FF Integration Fix)
- Re-read `task_plan.md` before edits and targeted the production-stage backbone FF activation path in workflow scripts.
- Fixed `run_sim_1rkl.sh` production wiring:
  - `prepare_hybrid_artifacts()` now calls `prepare_upside_backbone_reference`.
  - production `augment_production_rotamer_nodes` call now passes `${UPSIDE_BACKBONE_REF_UP}` (required 4th arg).
- Fixed `test_prod_run_sim_1rkl.sh` to match production backbone integration requirements:
  - Added Upside backbone reference config/validation (`UPSIDE_BACKBONE_REF_*`, rama/hbond/reference-state paths).
  - Added `prepare_upside_backbone_reference()` and invoked it before stage 7.0 preparation.
  - Updated production node injection to require/reference the generated backbone `.up` and to copy/remap canonical backbone nodes (`Distance3D`, `Angle`, `Dihedral_*`, `Spring_*`, `rama_*`, `infer_H_O`, `protein_hbond`, `hbond_energy`, `backbone_pairs`).
  - Updated production augment call to pass `${UPSIDE_BACKBONE_REF_UP}`.
- Validation:
  - `bash -n run_sim_1rkl.sh` passed.
  - `bash -n test_prod_run_sim_1rkl.sh` passed.
  - Smoke run: `PROD_70_NSTEPS=1 PROD_FRAME_STEPS=1 bash test_prod_run_sim_1rkl.sh` completed stage prep/run without missing-node failure.
  - HDF5 verification on `outputs/martini_test_1rkl_hybrid/checkpoints/1rkl.stage_7.0.up` confirmed:
    - required backbone nodes present in `/input/potential`.
    - remapped backbone atom IDs are within runtime AA carrier range (`reference_index_offset=4216`, `reference_index_count=233`).
    - stage atom count includes carriers (`n_atom=4449`).

## 2026-02-19 (Rotamer-State Read Path Audit)
- Re-read `task_plan.md` and audited production hybrid rotamer coupling path in C++.
- Verified runtime binding/scheduling in `../../src/martini.cpp`:
  - `register_hybrid_for_engine()` discovers/binds `rotamer` and placement nodes from live DerivEngine graph.
  - Adds explicit dependency edges so `martini_potential*` depends on placement + rotamer + coord nodes.
  - In production active mode, startup fails hard if rotamer/placement nodes are missing.
- Verified live probability pull each MARTINI force evaluation:
  - `MartiniPotential::compute_value()` calls `refresh_sc_row_probabilities_from_rotamer()` every call.
  - That function reads `rotamer_node->get_value_by_name("node_marginal")` and `("node_lookup")` (no static fallback).
  - If lookup/marginal mapping is incomplete or invalid, it returns false and production throws a fatal error.
- Verified rotamer node semantics in `../../src/rotamer.cpp`:
  - `compute_value(DerivMode|PotentialAndDerivMode)` solves marginals (`solve_for_marginals()`) each call.
  - `get_value_by_name("node_marginal")` returns current `cur_belief` marginals from that solve.
- Verified integration ordering behavior:
  - Default `--integrator v` path calls `DerivEngine::compute(DerivMode)` each sub-step.
  - Added dependencies ensure rotamer executes before martini potential in topological order for each sub-step.

## 2026-02-19 (Single-File Upside Backbone FF Integration)
- User requested removing separate backbone reference `.up` generation and integrating Upside backbone FF node generation directly into hybrid production-stage preparation.
- Updated architecture record in `task_plan.md` Revised Decisions:
  - production stage now generates backbone FF nodes in-place in the stage `.up` via `upside_config` writer functions (single-file workflow).
- Implemented in both workflow scripts:
  - `run_sim_1rkl.sh`
  - `test_prod_run_sim_1rkl.sh`
- `augment_production_rotamer_nodes` now:
  - loads `upside_config.py` dynamically from `${UPSIDE_HOME}/py/upside_config.py`;
  - creates backbone nodes directly in current stage file (`Distance3D`, `Angle`, `Dihedral_*`, `Spring_*`, `rama_coord`, `rama_map_pot`, `rama_map_pot_ref`, `infer_H_O`, `protein_hbond`, `hbond_energy`, `backbone_pairs`) using writer functions;
  - remaps generated protein-local backbone atom indices to runtime AA carrier indices in the same stage file.
- Removed secondary backbone-reference workflow usage from scripts:
  - no call to `prepare_upside_backbone_reference` in execution path;
  - removed corresponding variables/function definitions from both scripts to avoid hidden side path.
- Validation:
  - `bash -n run_sim_1rkl.sh` passed.
  - `bash -n test_prod_run_sim_1rkl.sh` passed.
  - Smoke run passed: `PROD_70_NSTEPS=1 PROD_FRAME_STEPS=1 bash test_prod_run_sim_1rkl.sh`.
  - Stage-7 logs confirm hybrid production active with rotamer/placement nodes, and stage runs without requiring a separate backbone reference `.up` file.

## 2026-02-19 (Input `.up` Visualization Script Fix)
- Re-read `task_plan.md` and debugged user-saved `visualize_up_input.py`.
- Reproduced immediate runtime failure:
  - `SyntaxError: unterminated string literal` caused by wrapped lines in f-strings and function calls.
- Patched `visualize_up_input.py`:
  - fixed broken wrapped lines (`residue_ids` assignment, charge summary print, scatter calls).
  - made input argument optional with default `inputs/1rkl.up`.
  - kept direct executable invocation path (`./visualize_up_input.py`) by using project venv shebang `#!../../.venv/bin/python`.
  - added writable matplotlib/fontconfig cache setup via `MPLCONFIGDIR` and `XDG_CACHE_HOME`.
  - auto-selects `Agg` backend when running headless (`DISPLAY` absent), so PNG export works non-interactively.
- Validation:
  - `./visualize_up_input.py --output /tmp/1rkl_input_default.png` passed and wrote image.
  - `./visualize_up_input.py inputs/1rkl.up --color-by type --point-size 4 --output /tmp/1rkl_input_type.png` passed and wrote image.
  - script remains executable: `-rwxr-xr-x`.

## 2026-02-19 (Visualization Script Save-Only Behavior)
- Updated `visualize_up_input.py` to always save figures to file instead of calling `plt.show()`.
- Added default output path behavior when `--output` is omitted:
  - saves to `<input_basename>_input_viz.png` in the same directory as the input `.up`.
- Validation:
  - `./visualize_up_input.py` now saves to `inputs/1rkl_input_viz.png`.
  - no non-interactive `FigureCanvasAgg` show warning is emitted.

## 2026-02-19 (Production Fallback Removal + Stability Debug/Fix)
- Re-read `task_plan.md`, audited `run_sim_1rkl.sh`/`test_prod_run_sim_1rkl.sh` and reproduced production instability:
  - `PROD_70_NSTEPS=1 PROD_FRAME_STEPS=1 ./test_prod_run_sim_1rkl.sh`
  - observed stage-7 jump (`martini_potential` from `~2.18e6` initial to `~4.80e7` at first MD step).
- Identified remaining production fallback path:
  - `/input/hybrid_sc_map/proj_target_indices` in stage-7 still pointed to MARTINI `BB` atoms.
  - runtime `project_sc_gradient_if_active()` still had BB fallback when explicit SC targets were absent.
- Patched workflow stage-file injection in both scripts:
  - `run_sim_1rkl.sh`
  - `test_prod_run_sim_1rkl.sh`
  - `inject_hybrid_mapping()` now rewrites SC projection targets/weights to runtime AA carrier indices (`N/CA/C/O`) using rewritten `hybrid_bb_map` component mappings.
  - adds SC mapping attrs:
    - `target_index_space=stage_runtime`
    - `target_projection=bb_component_carriers`
- Patched runtime strictness in `../../src/martini.cpp`:
  - removed SC->BB fallback in `project_sc_gradient_if_active()`.
  - added strict SC target validation:
    - each SC row must have explicit nonzero projection targets;
    - SC targets cannot be MARTINI `BB`/`SC` role atoms.
- Rebuilt C++:
  - `cmake --build ../../obj -j4` succeeded.
- Found additional root cause for immediate stage-7 instability:
  - stage handoff left injected AA carrier coordinates stale when copying positions from stage 6.6 to stage 7.0.
  - measured BB/carrier mismatch after handoff was ~`55.38 Å` RMS (`/tmp` probe).
- Patched `set_initial_position.py`:
  - added `refresh_hybrid_reference_carriers(...)` to realign injected carrier coordinates so each residue’s carrier COM matches current BB proxy position during handoff.
  - post-fix mismatch probe: ~`3e-7 Å` RMS.
- Post-fix validation:
  - `PROD_70_NSTEPS=1 PROD_FRAME_STEPS=1 ./test_prod_run_sim_1rkl.sh`
    - first-step jump removed (`martini_potential` remained at initial scale `~2.21e6`).
  - confirmed stage-7 SC projection targets are AA-only (`N/CA/C/O`), zero `BB`/`SC` targets.
- Continued debug for remaining multi-step instability:
  - with production timestep `0.020`, 20-step run still blew up later (`~step 12+`) even with SC rows disabled.
  - isolated root cause: timestep mismatch after injecting stiff Upside all-atom backbone nodes into production stage.
- Fixed production timestep defaults:
  - `run_sim_1rkl.sh`: `PROD_TIME_STEP` default changed `0.020 -> 0.002`.
  - `test_prod_run_sim_1rkl.sh`: `PROD_TIME_STEP` default changed `0.020 -> 0.002`.
- Stability validation with new defaults/fixes:
  - `bash -n run_sim_1rkl.sh` and `bash -n test_prod_run_sim_1rkl.sh` passed.
  - `PROD_70_NSTEPS=5 PROD_FRAME_STEPS=1 ./test_prod_run_sim_1rkl.sh` ran stably with `dt=0.002` and finite energies.
  - direct 20-step probes at `dt=0.002` (with and without SC rows) remained stable and finite (no blow-up/NaN).

## 2026-02-19 (Preparation Alignment: AA Backbone COM -> MARTINI BB)
- Implemented preparation-time RMSD alignment in `prepare_hybrid_system.py`:
  - Added `compute_backbone_com_alignment(...)` to fit a rigid transform (Kabsch) from all-atom residue backbone COM points (`N/CA/C/O`) to MARTINI `BB` positions.
  - `collect_bb_map(...)` now applies this transform when writing `hybrid_bb_map/reference_atom_coords`.
  - Added alignment diagnostics in preparation output and `bb_comment` (`align_rmsd=...`).
- Validation:
  - `python3 -m py_compile prepare_hybrid_system.py` passed.
  - Ran generation:
    - `python3 prepare_hybrid_system.py --protein-pdb pdb/1rkl.pdb --protein-cg-pdb pdb/1rkl_hybrid.MARTINI.pdb --protein-itp pdb/1rkl_hybrid_proa.itp --bilayer-pdb pdb/bilayer.MARTINI.pdb --output-dir /tmp/hybrid_align_check`
  - Runtime output showed:
    - `Backbone COM RMSD alignment (AA->MARTINI BB): n=31, rmsd=0.1030 Å`
  - HDF5 check confirmed aligned coordinates are stored in `/input/hybrid_bb_map/reference_atom_coords` and comments include `align_rmsd=0.1030`.

## 2026-02-19 (Rotamer Projection + Probabilistic SC-BB Weighting Audit/Fix)
- Re-read `task_plan.md`, audited runtime hybrid coupling and stage-7 artifacts for:
  - rotamer projection anchor frame (all-atom backbone carriers),
  - live probabilistic weighting path.
- Verified stage-7 data integrity on `outputs/martini_test_1rkl_hybrid/checkpoints/1rkl.stage_7.0.up`:
  - `/input/potential/affine_alignment/atoms` uses `N/CA/C` carrier roles only.
  - `/input/hybrid_sc_map/proj_target_indices` nonzero targets are `N/CA/C/O` only.
  - required nodes present: `rotamer`, `placement_fixed_point_vector_only`, `placement_fixed_scalar`, `affine_alignment`.
- Verified live rotamer marginal source:
  - `../../src/rotamer.cpp` exports `node_marginal` from `cur_belief` after `solve_for_marginals()`.
  - `../../src/martini.cpp` refreshes SC row probabilities from `rotamer` node each MARTINI compute call.
- Found and fixed a gap:
  - before fix, probabilistic path covered only protein-vs-nonprotein SC interactions; allowed same-residue protein `SC-BB` pairs took deterministic MARTINI pair path.
  - patched `../../src/martini.cpp` to route allowed same-residue protein `SC-BB` pairs through the same probabilistic per-rotamer weighting path.
  - added frame-correct env-gradient back-rotation for protein BB partners when coupling alignment is enabled.
- Rebuilt binary:
  - `cmake --build ../../obj -j 8` (success).
- Validation:
  - short production smoke passed: `PROD_70_NSTEPS=5 PROD_FRAME_STEPS=5 ./test_prod_run_sim_1rkl.sh`.
  - stage-7 pair-table audit confirms patched path is active for real data:
    - same-residue protein `SC-BB` pairs eligible for probabilistic weighting: `6`.

## 2026-02-19 (Production Potential Fluctuation Root-Cause Debug + Fix)
- Reproduced user-reported production-stage pathological MARTINI potential swings on controlled 200-step probes from the same stage-7 initial state:
  - full SC probabilistic coupling (before fix): `martini_potential` `2.229e6 -> 2.443e7 (step 50) -> 1.215e5 (step 100) -> 1.509e4 (step 150)`.
- Ruled out non-protein hard-sphere branch as primary source:
  - toggled `production_nonprotein_hard_sphere` (`1` vs `0`) from identical initial checkpoint; both showed the same early large spikes.
- Isolated source to SC probabilistic energy aggregation:
  - disabling SC map rows (`n_sc=0`) removed the pathological regime immediately (`martini_potential` stayed `~2e3-3e3` in 200 steps).
  - this showed the instability comes from probabilistic SC interaction accumulation, not the base MARTINI pair engine.
- Root cause identified:
  - SC probabilistic coupling used linear expectation over rotamer-state energies (`sum_r p_r E_r`) for each proxy-env pair.
  - low-probability high-clash rotamers therefore contributed disproportionately large positive energies, destabilizing long production runs.
- Implemented fix in `../../src/martini.cpp`:
  - replaced linear per-edge averaging with Boltzmann-consistent free-energy mixing:
    - potential contribution: `-log(sum_r p_r * exp(-E_r))`
    - gradient weights: posterior `p_r * exp(-E_r)` normalized per edge.
  - added finite-math-safe stabilization:
    - shifted energies by per-edge minimum,
    - clamped exponent argument range,
    - floor on partition sum and finite checks with prior-weighted fallback only for degenerate numeric cases.
- Rebuilt binary:
  - `cmake --build ../../obj -j 8` passed.
- Validation (post-fix):
  - 200-step probe from same initial stage-7 state:
    - `martini_potential` `4.687e3 (step 0)`, `3.812e3 (50)`, `2.912e3 (100)`, `4.103e3 (150)`, all finite.
  - 600-step probe from fresh stage-7 state:
    - stable through reported trouble window (`step 500`: `Rg=13.0 Å`, `martini_potential=1.147e4`, no blow-up/NaN).

## 2026-02-20 (Production Instability Follow-up: SC Relaxation Sign Fix)
- Re-read `task_plan.md` and investigated user-reported stage-7 instability (`martini_potential` spike near step 1000 with rapid hbonds loss).
- Confirmed latest binary now includes hard-sphere cap controls in startup log:
  - `hs_force_cap=100.000`, `hs_pot_cap=5000.000`.
- Found and fixed a sign bug in `../../src/martini.cpp` SC inner relaxation loop:
  - previous logic advanced SC proxy positions with `+dE/dx` (gradient ascent),
  - patched to use force direction `-dE/dx` for the relaxation update.
- Rebuilt C++:
  - `cmake --build ../../obj -j4` (success).
- Diagnostics:
  - exact production settings (`sc_env_relax_steps=200`) now start from lower MARTINI energy on same input (`4686.84 -> 3436.05` at step 0).
  - short fast diagnostic with same code path but `sc_env_relax_steps=20` and frequent reporting completed stably for 400 steps:
    - step 50: `martini_potential=2980.32`
    - step 100: `2444.70`
    - step 200: `3084.59`
    - step 300: `4886.27`
    - step 350: `4496.45`
    - no NaN/overflow observed in this probe.
- Remaining work:
  - long-horizon validation under exact production settings (`sc_env_relax_steps=200`, 5k steps) is still pending/running due runtime cost.

## 2026-02-20 (Production Start Protein Placement Fix)
- Investigated user-reported production-start artifact: protein appears at box corner and is visually PBC-split at stage 7 start.
- Root cause identified in stage handoff:
  - `set_initial_position.py` copied last frame positions from stage 6.6 to stage 7.0 without any production recentering.
  - In checkpoint data, protein coordinates were near origin while box center was near `(Lx/2, Ly/2, Lz/2)`, so wrapped visualization placed protein near boundaries.
- Implemented fix in `set_initial_position.py`:
  - added `recenter_protein_for_production(...)` executed only when stage label is `production`;
  - uses `/input/hybrid_bb_map/bb_atom_index` as protein reference set;
  - computes protein BB center in periodic space via circular mean per axis;
  - shifts all coordinates so BB center maps to box center, then wraps coordinates into `[0, L)`.
- Validation:
  - syntax check: `python3 -m py_compile set_initial_position.py` passed.
  - handoff test:
    - `python set_initial_position.py outputs/.../1rkl.stage_6.6.up /tmp/1rkl.stage_7.0.recenter_test.up`
    - script log reported recenter action and updated box.
  - post-check on `/tmp/1rkl.stage_7.0.recenter_test.up`:
    - BB circular center matches target box center to numerical precision.
    - BB coordinates are compact around box center rather than boundary-origin.
