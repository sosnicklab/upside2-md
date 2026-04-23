# Findings

## 2026-04-23 (`run_sim_1afo_outlipid.sh` Isolation Verification)
- The Slurm/bootstrap fix in `run_sim_1afo.sh` does not alter the intended continuation behavior of `run_sim_1afo_outlipid.sh`.
- Verified behavior after the AABB wrapper fix:
  - `run_sim_1afo_outlipid.sh` still exports `DISABLE_1AFO_AABB_AUTO_CONTINUE=1` before delegating,
  - when only an AABB prior stage-7 artifact exists, the outlipid wrapper does not resume,
  - when an outlipid prior stage-7 artifact exists, the outlipid wrapper resumes from that file and uses the outlipid continue directory,
  - under a simulated `sbatch` context with a wrong inherited `UPSIDE_HOME`, the outlipid wrapper still delegates with the correct resolved project root.
- Consequence:
  - the AABB Slurm fix is isolated to the AABB wrapper behavior and does not regress the outlipid wrapper.

## 2026-04-23 (`run_sim_1afo.sh` Slurm Bootstrap Failure Root Cause)
- The `sbatch run_sim_1afo.sh` failure with
  - `ERROR: UPSIDE executable not found: /home/yinhanw/project/yinhan/upside2-md/obj/upside`
  was caused by the wrapper missing the Slurm-safe bootstrap that already existed in the outlipid wrappers.
- The previous `run_sim_1afo.sh`:
  - did not search for `run_sim_1rkl.sh` using `SLURM_SUBMIT_DIR`,
  - did not load the cluster modules / activate the repo `.venv`,
  - did not force `UPSIDE_SKIP_SOURCE_SH=1`,
  - did not force `UPSIDE_HOME` to the resolved project root before delegating.
- Consequence:
  - `sbatch` could inherit a stale or unrelated `UPSIDE_HOME`,
  - the delegated base workflow then looked for `obj/upside` under that wrong root and aborted before any simulation work started.
- Correct fix:
  - make `run_sim_1afo.sh` use the same Slurm-safe bootstrap pattern as `run_sim_1rkl_outlipid.sh` / `run_sim_1afo_outlipid.sh`.
- Verification after the fix:
  - isolated Slurm-style harness with `UPSIDE_HOME=/wrong/root` still delegated with:
    - `UPSIDE_HOME=<resolved project root>`
    - `UPSIDE_SKIP_SOURCE_SH=1`
    - `PYTHONPATH=<resolved project root>/py`
    - `PATH` containing `<resolved project root>/obj`

## 2026-04-23 (1AFO Sidechain Injection Failure Root Cause)
- The `1afo_outlipid` failure in `inject_stage7_sc_table_nodes(...)` was caused by a multi-chain residue-ID collapse in the shared hybrid metadata path, not by the `1afo` wrapper scripts themselves.
- For `example/16.MARTINI/pdb/1AFO.pdb`:
  - `extract_backbone_sequence(...)` correctly returns `72` residues (`36` in chain `A`, `36` in chain `B`),
  - `collect_aa_backbone_map(...)` previously stored raw PDB `resseq` as `hybrid_bb_map/bb_residue_index`,
  - because both chains use residue numbers `66..101`, the metadata contained only `36` unique backbone residue IDs.
- Consequence:
  - sidechain injection saw `72` sequence entries but only `36` affine residues from `build_affine_atoms(...)`,
  - which triggered:
    - `ValueError: Missing or inconsistent /input/sequence for AA-backbone sidechain injection: expected 36 residues`
- Correct fix:
  - assign a unique residue-order `bb_residue_index` per backbone residue during hybrid metadata generation,
  - do not use raw `resseq` alone as the unique hybrid residue key for multi-chain proteins.
- Verification after the fix:
  - helper output for `1AFO`: `72` sequence residues and `72` unique `bb_residue_index` values,
  - written metadata file `/tmp/1afo_test_backbone_metadata.h5`: `sequence len = 72`, `unique bb_residue_index = 72`,
  - reduced local `run_sim_1afo_outlipid.sh` smoke advanced past the prior `inject-stage7-sc` failure site.

## 2026-04-22 (`run_sim_1rkl_outlipid.sh` Validity Audit After Shared Workflow Updates)
- The current `run_sim_1rkl_outlipid.sh` does inherit the recent shared workflow/engine updates because it now delegates to `example/16.MARTINI/run_sim_1rkl.sh`.
- That inheritance is incomplete for the out-of-bilayer intent:
  - the wrapper exports `PROTEIN_PLACEMENT_MODE`, `PROTEIN_ORIENTATION_MODE`, and `PROTEIN_SURFACE_GAP`,
  - the current prep path in `py/martini_prepare_system.py::prepare_mixed_structure(...)` and `py/martini_prepare_system_lib.py` does not read those settings,
  - the wrapper also exports `PROTEIN_LIPID_MIN_GAP` and `PROTEIN_LIPID_CUTOFF_MAX`, but the current base workflow redefines those as local readonly constants (`4.5` and `8.0`) instead of reading env overrides.
- A direct prep check with the wrapper's current box/padding/cutoff values confirmed the functional regression:
  - generated file: `example/16.MARTINI/outputs/outlipid_validity_check/prep/1rkl_outlipid_check.MARTINI.pdb`
  - protein `z_min = 55.143 Å`
  - upper-leaflet `PO4 z_max = 91.635 Å`
  - clearance above upper leaflet = `-36.492 Å`
- Consequence:
  - the current script is shell-valid and picks up the new rigid-body preproduction behavior from the base workflow,
  - but it is not valid as an "initially outside of bilayer" workflow because the generated starting geometry remains bilayer-embedded.

## 2026-04-21 (1RKL AABB Stage-7 Oxygen Failure Root Cause)
- The reported bilayer "shifting through space" in `outputs/martini_test_1rkl_aabb` is not the dominant runtime bug.
- Direct artifact measurements showed:
  - preproduction `PO4` motion relative to the fixed protein stays below about `0.15 Å` in XY and `0 Å` in Z,
  - stage-7 absolute `PO4` drift is small (`~0.31 Å` XY),
  - the large apparent stage-7 motion comes from the protein destabilizing relative to a mostly stable bilayer.
- Root cause:
  - `py/martini_prepare_system_lib.py::inject_backbone_nodes(...)` writes the native Upside 3-site backbone machinery for `N/CA/C`,
  - explicit runtime `O` atoms are present in the AA-backbone workflow but are not included in those generated backbone bond/angle nodes,
  - stages `6.0` through `6.6` hid this because `fix_rigid` held the whole AA backbone fixed,
  - stage `7.0` releases the backbone, so the unrestrained `O` atoms fly off immediately.
- Correct fix for the current workflow:
  - keep the existing 3-site Upside backbone nodes,
  - append explicit runtime oxygen geometry terms in the stage file:
    - `C-O` bond,
    - `CA-C-O` angle,
    - `O-C-N(next)` angle for nonterminal residues.
- Verification on a regenerated production handoff confirmed the fix:
  - original `stage_7.0.up`: `C-O` max grew to `102.05 Å`,
  - fixed stage-7 smoke run over `1000` MD steps: `C-O` stayed within `0.96–1.58 Å`.

## 2026-04-06 (Protein-Backbone RMSD Alignment Audit)
- There is no generic Upside MD step that RMSD-aligns and overwrites protein backbone coordinates during integration.
- The relevant feature is MARTINI-hybrid specific:
  - runtime state includes `integration_rmsd_align_enable` in `src/martini.cpp`,
  - `example/16.MARTINI/run_sim_1rkl.sh` writes that attr into production hybrid control,
  - the workflow default is `INTEGRATION_RMSD_ALIGN_ENABLE=1`.
- `martini_hybrid::align_active_protein_coordinates(...)` is called from `src/deriv_engine.cpp` before derivative evaluations inside the integrator, including multiple sub-stages per MD step.
- Current behavior of that function is bookkeeping only:
  - it refreshes reference backbone coordinates,
  - computes a Kabsch/Horn transform against the previous frame/reference,
  - stores `prev_bb_pos_rmsd`,
  - increments `integration_align_step`,
  - optionally prints debug info.
- It does not apply the transform to `pos->output` or `mom`, so it does not rotate/translate the actual integrated trajectory.
- The existing `affine_alignment` node elsewhere in Upside is a residue-local frame construction node, not a whole-protein per-step RMSD trajectory alignment.

## 2026-04-06 (Would Removing Whole-Protein Rigid-Body Momentum Change Results?)
- Yes, it would change the simulation dynamics.
- Removing total linear COM momentum is common and typically does not materially change translationally invariant configurational observables, but it still changes the lab-frame trajectory and momentum distribution.
- Removing whole-protein angular momentum is more invasive:
  - it suppresses rigid-body rotational diffusion,
  - it changes the angular-momentum distribution generated by the thermostat,
  - if applied only to the protein in a mixed system, it also changes protein-environment momentum exchange.
- Therefore angular-momentum removal is not a "no-change" fix for Upside results.
- If the observed rolling is caused by a real force/torque imbalance, the correct fix is to remove that imbalance.
- If the rolling is just thermostat-driven rigid-body motion, the cleanest non-invasive path is to leave dynamics alone and compute diagnostics or coupling quantities in a co-moving/body frame rather than altering simulation coordinates or momenta.

## 2026-04-06 (Workflow Rollback: Disable Hybrid RMSD-Align)
- Flipping only `INTEGRATION_RMSD_ALIGN_ENABLE` in `run_sim_1rkl.sh` was insufficient.
- Root cause:
  - early hybrid stages inherit `/input/hybrid_control` from the generated mapping/prep artifacts,
  - when `integration_rmsd_align_enable` is absent, `src/martini.cpp` reads it with default `1`.
- Therefore the correct workflow-level rollback point is the hybrid-control metadata writer in `py/martini_prepare_system_lib.py`, not only the later production-stage override in `run_sim_1rkl.sh`.
- After writing `integration_rmsd_align_enable=0` at metadata creation time, a reduced end-to-end workflow run completed through stage `7.0` with runtime banners showing `integration_rmsd_align=0`.

## 2026-04-06 (Full Feature Removal: Hybrid RMSD-Align)
- The bookkeeping path was fully removable:
  - no other runtime feature depended on `align_active_protein_coordinates(...)`,
  - the Kabsch/Horn helpers in `src/martini.cpp` were dead once that function was removed,
  - the integrator ran normally after deleting all five call sites in `src/deriv_engine.cpp`.
- The workflow also no longer needs any RMSD-align configuration surface:
  - `INTEGRATION_RMSD_ALIGN_ENABLE` in `run_sim_1rkl.sh` was deleted,
  - `integration_rmsd_align_enable`, `coupling_align_debug`, and `coupling_align_interval` are no longer written into hybrid-control metadata.
- Reduced end-to-end verification after rebuild completed through stage `7.0` with runtime banners that no longer mention any RMSD-align field at all.

## 2026-04-09 (Thermostat Default Audit For Workflow 16)
- The standard Upside examples do not pass `--thermostat-timescale`:
  - `example/01.GettingStarted/0.run.py`,
  - `example/08.MembraneSimulation/0.normal.run.py`,
  - `example/12.MultistepIntegrator/0.run.py`.
- Therefore they inherit the engine default thermostat timescale `tau = 5.0` from `src/main.cpp`.
- Workflow `16` had been the only example in the main example set that explicitly overrode this to `THERMOSTAT_TIMESCALE=4.0` in `run_sim_1rkl.sh`, with no local documentation explaining why.
- To keep the damping convention consistent with the examples that underlie the user’s `40 ps` per-step calibration, workflow `16` should default to `tau = 5.0` as well.
- Separate caveat:
  - workflow `16` still uses `EQ_TIME_STEP=0.010` and `PROD_TIME_STEP=0.002`,
  - the standard examples rely on the engine default `dt = 0.009`,
  - so matching the thermostat default does not, by itself, make the borrowed `40 ps/step` calibration valid for hybrid stage `7.0`.

## 2026-04-14 (Hybrid Packing Analysis Surface)
- The current hybrid workflow already carries production-stage NPT controls, but leaves them off by default:
  - `PROD_70_NPT_ENABLE=0` in `example/16.MARTINI/run_sim_1rkl.sh`.
- The stage-7 output schema is sufficient to measure bilayer packing directly from local outputs:
  - `output/pos`,
  - `output/time`,
  - `output/box` when NPT writes it,
  - otherwise static `input/potential/martini_potential.{x_len,y_len,z_len}` attrs.
- The current hybrid stage files also contain enough lipid topology to compute packing observables without external topology files:
  - lipid molecules can be identified from `PO4` atom names plus `molecule_ids`,
  - current DOPC tail bead names are present (`GL1/GL2`, `C1A/C2A/D3A/C4A/C5A`, `C1B/C2B/D3B/C4B/C5B`).
- Consequence:
  - the missing piece for a packing-first calibration loop was not runtime support,
  - it was a direct analysis path for APL / thickness / order on completed hybrid outputs.

## 2026-04-14 (Current Hybrid Packing Baseline)
- The new packing analyzer ran successfully on:
  - `example/16.MARTINI/outputs/martini_test_1rkl_hybrid/checkpoints/1rkl.stage_7.0.up`
- Measured baseline after `20%` burn-in:
  - `n_lipid_molecules = 280`
  - `APL = 89.01 A^2`
  - `PO4 thickness = 32.57 A`
  - `tail order parameter = 0.238`
- That output was produced without dynamic `output/box` frames.
  - The analyzer therefore fell back to the static `input/potential/martini_potential.{x_len,y_len,z_len}` box attrs and reported zero area variance, which is the expected NVT behavior for this artifact.
- Consequence:
  - the repository can now quantify current hybrid membrane packing directly,
  - the next scientific decision still needs an external dry-MARTINI packing target for comparison rather than more blind interface-scale tuning.

## 2026-04-19 (1RKL Mixed-System Placement Audit For Out-of-Bilayer Start)
- The current mixed-system prep path in `py/martini_prepare_system.py::prepare_mixed_structure(...)` does not have an explicit outside-of-bilayer mode.
- Current behavior:
  - translates the protein center to the bilayer center,
  - sizes the XY target from the translated protein span,
  - removes overlapping lipids,
  - recenters the bilayer in `z` when defining the box.
- Consequence:
  - an out-of-bilayer starting geometry cannot be implemented reliably in a wrapper script alone,
  - the correct change point is the mixed-system prep step that writes the packed runtime PDB and hybrid mapping inputs.
- The `1rkl` input PDB is already membrane-oriented metadata-wise (`REMARK 1/2 of bilayer thickness: 12.4`), but a reviewable "protein on top of bilayer" start still needs an explicit orientation rule rather than only a `z` translation.

## 2026-04-19 (Out-of-Bilayer 1RKL Review Geometry)
- The new outside-start placement that best matches the request is:
  - `protein_placement_mode = outside-top`
  - `protein_orientation_mode = lay-flat`
  - `protein_surface_gap = 6.0 Å`
- A direct stage-0 preparation run with those settings produced:
  - box = `177.472 x 177.472 x 160.621 Å`
  - protein span after flat orientation = `43.177 x 16.204 x 11.812 Å`
  - top-of-bilayer clearance = exactly `6.0 Å`
  - removed lipid residues = `0`
- Consequence:
  - the protein is above the upper leaflet with explicit surface separation,
  - the larger box comes primarily from the widened XY target plus increased `z` padding,
  - the start geometry does not rely on carving a hole into the bilayer to make room.

## 2026-04-19 (Stage-7 Continuation Path For Workflow 16)
- The existing helper `set-initial-position` already supports production-stage continuation because it can copy the last `/output/pos` frame from one `.up` file into another `.up` file's `/input/pos` while preserving the last box.
- `martini_extract_vtf.py` can extract directly from a resumed `stage_7.0.up` without requiring regenerated stage-0 prep artifacts.
- Consequence:
  - the clean restart design is to skip stage `0` through `6.6`,
  - validate the previous `stage_7.0.up` is still a hybrid production file,
  - write continuation output to a new `stage_7.0.continue.up` by default,
  - then run production MD and VTF extraction from that resumed file.

## 2026-04-19 (User Correction: Slurm Walltime Hard Limit)
- This cluster's hard Slurm walltime limit for the workflow is `36:00:00`, not `48:00:00`.
- Lesson:
  - do not invent or round up Slurm walltimes in wrapper scripts;
  - when the user provides a hard cluster limit, propagate that exact value into `#SBATCH --time` and any related guidance.

## 2026-04-19 (User Correction: Slurm Wrapper Path Resolution)
- `sbatch` may execute a spool copy of the submitted script under `/var/spool/slurm/...`, so `BASH_SOURCE[0]` is not a reliable way for a wrapper to find sibling workflow files.
- Lesson:
  - do not resolve repo-local helper scripts in an `sbatch` wrapper from the wrapper's own runtime path alone;
  - prefer `SLURM_SUBMIT_DIR`, the job working directory, or an explicit override env var such as `BASE_WORKFLOW_SCRIPT`.

## 2026-04-19 (User Correction: Proper Slurm Environment For This Project)
- The proper cluster pattern for this project must be documented directly in the main repo guidance, not inferred from side directories that may be removed.
- The self-contained Slurm environment contract is:
  - source `/etc/profile.d/modules.sh` when present,
  - `module load python/3.11.9`,
  - `module load cmake`,
  - `module load openmpi`,
  - `module load hdf5/1.14.3` or the configured override,
  - source `PROJECT_ROOT/.venv/bin/activate` when present,
  - export `UPSIDE_HOME`, `PATH`, and `PYTHONPATH` from the repo root.
- Lesson:
  - for Slurm jobs, do not depend on the repo's Apple-Silicon `source.sh` bootstrap to establish the runtime environment;
  - Slurm wrappers should establish the cluster environment themselves and set `UPSIDE_SKIP_SOURCE_SH=1` before invoking lower-level workflow scripts;
  - `AGENTS.md` must describe the Slurm setup directly and must not tell future instances to look at removable auxiliary folders for core environment rules.

## 2026-04-20 (User Correction: VTF Output Must Be Per Trajectory Segment)
- A continued `.up` file may contain multiple trajectory groups:
  - archived `output_previous_*` groups from earlier segments,
  - the current `output` group for the latest segment.
- The correct VTF export behavior is one file per trajectory segment, not a flattened output that combines all archived pieces into one `.vtf`.
- Lesson:
  - when a restartable HDF5 stage file can archive prior outputs, inspect the internal group layout before assuming one stage file implies one trajectory segment;
  - keep the single-segment filename stable, but split multi-segment extraction into deterministic `segment_<n>` files.

## 2026-04-20 (User Correction: Workflow Seeds Must Not Be Fixed Literals)
- `example/16.MARTINI/run_sim_1rkl.sh` still carried fixed default literals for:
  - `PREP_SEED=2026`
  - `SEED=7090685331`
- The correct default behavior for this workflow is:
  - generate fresh seeds per run when the user does not provide them,
  - keep explicit env-var overrides available for reproducible reruns.
- Lesson:
  - when a workflow is intended for independent reruns, audit seed defaults directly instead of assuming a variable named `SEED` is already randomized;
  - avoid fixed literal seed defaults unless the user explicitly asked for deterministic behavior.

## 2026-04-20 (Audit: 1RKL Workflow Integrator Path)
- `example/16.MARTINI/run_sim_1rkl.sh` uses `--integrator v` for both minimization and MD stages rather than a dedicated MARTINI-specific CLI mode.
- `src/main.cpp` always calls `martini_masses::load_masses_for_engine(...)` after loading the HDF5 config, so MARTINI mass data is available whenever the stage file contains `/input/mass`.
- The generated 1RKL stage files do contain `/input/mass`; inspection of existing workflow outputs confirms that both stage `6.0` and stage `7.0` checkpoints include that dataset.
- Therefore the workflow does use a mass-aware integrator path at runtime, but the active implementation is the inline `martini_masses::has_masses(...)` branches in `src/deriv_engine.cpp`.
- The helper `martini_masses::martini_integration_stage(...)` in `src/martini.cpp` is not called anywhere and is effectively unused by this workflow.

## 2026-04-20 (Cleanup: Remove Dead `martini_integration_stage`)
- The only occurrence of `martini_integration_stage(...)` was its unused definition in `src/martini.cpp`.
- Removing that helper does not affect the active MARTINI mass-aware path, because runtime integration still uses:
  - `martini_masses::load_masses_for_engine(...)` in `src/main.cpp`,
  - `martini_masses::has_masses(...)` and `martini_masses::get_mass(...)` inside `src/deriv_engine.cpp`.
- A fresh out-of-tree build from `src/CMakeLists.txt` completed successfully after the removal.

## 2026-04-21 (User Correction: Remove Retired Workflow-16 Code Completely)
- For the `1rkl` AA-backbone workflow replacement, it is not sufficient to block or ignore the old dry-MARTINI protein presentation path.
- The expected implementation standard is:
  - delete the retired workflow-16 CG/hybrid branches,
  - delete their dead control surface,
  - keep only code that remains active in the new design.
- Lesson:
  - when the user asks to retire a workflow path, remove the superseded implementation rather than leaving compatibility branches unless there is a concrete active caller that still depends on them;
  - if shared metadata names remain in use, that is fine, but the retired behavior behind them must still be deleted.

## 2026-04-21 (User Correction: Prevent Preproduction Bilayer Z Drift)
- The AA-backbone workflow cannot rely on soft preproduction behavior alone if bilayer `z` drift destabilizes stage `7.0`.
- For workflow `16`, the intended preproduction behavior is:
  - protein held rigid in absolute space,
  - bilayer `PO4` headgroups held at fixed `z` coordinates through the preproduction stages before production.
- Lesson:
  - when a user reports stage-to-stage geometric drift, verify whether the current mechanism is a hard hold or only a soft restraint;
  - if the request is for fixed coordinates, do not substitute a weak spring ramp when the engine can support an explicit fixed-axis hold.

## 2026-04-21 (User Correction: Preserve Committed Initial Bilayer Placement)
- For the `1rkl` prep path, it is not sufficient that the new code produces a reasonable bilayer placement.
- The initial bilayer structure and its `z` placement must remain identical to the last committed preparation behavior unless the user explicitly asks to change it.
- Lesson:
  - when replacing a representation inside an existing prep workflow, compare generated prep artifacts against the committed baseline instead of assuming equivalent centering code preserves the same structure;
  - if the user asks for identity with the committed baseline, verify the bilayer coordinates directly.

## 2026-04-21 (Workflow-16 Bilayer Asset Location And Provenance)
- The workflow-16 example should carry its default DOPC bilayer PDB locally under `example/16.MARTINI/pdb/DOPC.pdb` rather than pointing at the shared `parameters/dryMARTINI/` copy.
- The example README should state that the bilayer structure should be generated from CHARMM-GUI.

## 2026-04-21 (User Correction: README Must Not Prescribe `DOPC.pdb`)
- The README wording for workflow `16` must describe the bilayer input generically:
  - use the lipid bilayer structure downloaded from CHARMM-GUI,
  - do not imply the workflow requires the filename `DOPC.pdb`,
  - do not imply the bilayer identity must be DOPC unless the user explicitly says so.
- Lesson:
  - when documenting workflow inputs, distinguish between the current local default path used by scripts and the broader user requirement described in the README;
  - do not turn an implementation detail like a current filename into a workflow requirement unless the user explicitly asks for that contract.

## 2026-04-22 (User Correction: Keep Stage-Length Controls Public)
- For the `1rkl` workflow cleanup, the number of iterations / MD steps in each stage must remain part of the user-facing interface.
- It is acceptable to remove or hardcode internal prep/runtime/testing knobs, but not the per-stage run-length controls.
- Lesson:
  - when shrinking a workflow control surface, distinguish between true internal plumbing knobs and the run-length controls the user actively uses for experiment design;
  - do not collapse stage-specific iteration / step-count controls into fixed constants unless the user explicitly asks for that.

## 2026-04-22 (AA-Backbone Stage-7 Drift Diagnosis)
- Re-analysis of the current local artifacts shows the reported production-stage bilayer `z` shift is not a bulk bilayer translation bug.
- In `example/16.MARTINI/outputs/martini_test_1rkl_aabb/checkpoints/1rkl.stage_7.0.up`:
  - bilayer/environment/total `z` drift stays around `0.2–0.4 Å`,
  - protein center `z` drift reaches about `7 Å`,
  - the `6.6 -> 7.0` coordinate handoff is exact.
- The existing AA-backbone preproduction path holds the protein by absolute pinning:
  - `set_preproduction_spatial_holds()` writes `fix_rigid/atom_indices`,
  - `src/martini.cpp::apply_fix_rigid_*` zeroes the selected atoms' forces and momenta directly.
- User correction on intended semantics:
  - during minimization and preproduction the protein should remain rigid internally,
  - but it must not be locked in space,
  - rigid-body translation and rotation should both be allowed.
- Consequence:
  - this is not solvable by only rearranging stage timing or recentering flags;
  - it requires extending the `fix_rigid` runtime path with a rigid-body mode distinct from the existing absolute pin.
- Lesson:
  - distinguish "rigid internally" from "fixed in absolute coordinates";
  - when a user asks for the former, do not reuse an absolute freeze primitive without checking whether whole-body motion is supposed to remain allowed.

## 2026-04-22 (Rigid-Body Preproduction Hold Verification)
- Extending `fix_rigid` with `mode = rigid_body` was sufficient to keep the AA-backbone protein internally rigid during preproduction without using the old absolute pin path.
- Verified behavior on `example/16.MARTINI/example/16.MARTINI/outputs/martini_test_1rkl_rigidbody_smoke`:
  - `stage_6.2.up` and `stage_6.6.up` write `fix_rigid.mode = rigid_body`,
  - protein internal span changes stay at the `1e-6 Å` level,
  - `PO4` `z` holds remain exact,
  - the protein is no longer mathematically forced to zero absolute displacement; the short smoke shows small but nonzero COM motion.
- A direct 1000-step production continuation from the prepared stage-7 handoff stayed well behaved:
  - bilayer `z` drift max = `0.066 Å`,
  - protein `z` drift max = `0.052 Å`,
  - bilayer-minus-protein relative `z` change max = `0.118 Å`.
- Consequence:
  - the new rigid-body preproduction mode does not reintroduce bilayer bulk drift,
  - and it reduces the early production protein/bilayer relative drift compared with the old AABB artifact over the same 20-frame window.
