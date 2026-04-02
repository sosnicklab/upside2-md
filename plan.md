
## 2026-04-01 Stage-7 SC Table Integration

### Project Goal
- Consume the completed external `/Users/yinhan/Documents/SC-training` training results by building a native-unit `martini.h5`.
- Replace the legacy stage-7 sidechain back-mapping / rotamer production path in `example/16.MARTINI/run_sim_1rkl.sh` with a direct dry-MARTINI sidechain table interaction.
- Keep the new sidechain force field inactive for stages `6.0` through `6.6` and active only in production stage `7.0`.

### Architecture & Key Decisions
- `martini.h5` will be built from the assembled `sc_table.json` output and stored as a native-unit library (`nm`, `kJ/mol`) with a simple table schema:
  - residue order,
  - target dry-MARTINI type order,
  - uniform radial grid,
  - radial energy table.
- Simulation-side conversion into Upside units remains parameterized:
  - the stage-7 runtime node will read explicit `energy_conversion_kj_per_eup` and `length_conversion_angstrom_per_nm` attrs;
  - no numeric unit conversion will be baked into training outputs.
- Revised decision after user correction on 2026-04-02:
  - restore the reduced simulation-mass path (`ff_mass / 12`) in preparation and stage-file augmentation;
  - keep the weighted BB mapping revert in place;
  - retain the scientific-notation logging improvement only as a readability fix.
- Revised decision after user correction:
  - the stage-7 dry-MARTINI coupling must be integrated directly with Upside carriers, not evaluated on parallel protein MARTINI proxy coordinates and projected afterward.
- Sidechain/environment coupling will remain a separate stage-7 potential node, but it must be directly Upside-carried:
  - use deterministic Upside sidechain-bead coordinates from `affine_alignment -> placement_fixed_point_only_CB`;
  - act only between selected protein residues and non-protein dry-MARTINI particles;
  - add equal-and-opposite force to the dry particle and the Upside sidechain-bead coordinate so the backbone receives the feedback through the coord-node derivative chain;
  - all legacy protein MARTINI `SC` proxy nonbonded interactions in `martini_potential` must be disabled in stage 7.
- Backbone/environment coupling must also be directly Upside-carried:
  - treat the Upside `CA` carrier as the dry-MARTINI protein `BB` interaction site;
  - keep dry-MARTINI `BB` to surrounding non-protein dry particles active, but apply the resulting force directly to `CA` and the dry particle rather than to a proxy `BB` bead followed by projection;
  - revised after regression report: the old weighted `BB -> N/CA/C/O` active map has been restored for now; the attempted CA-only active BB remap has been reverted pending a deeper force-path audit;
  - revised after the 2026-04-02 production-drift audit:
    - direct `BB -> CA` and direct `CB` stage-7 forces must still honor the existing startup protein-feedback ramp (`sc_env_backbone_hold_steps`) instead of feeding full force back onto the protein from step `0`;
    - stage-7 must not evaluate any protein-internal MARTINI proxy-proxy terms, because those legacy bonded/nonbonded terms only add bookkeeping energy once the protein is carried directly by Upside;
  - do not evaluate protein internal SC-backbone dry-MARTINI interactions, because those are already handled on the Upside side.
- The workflow must remove the legacy production SC path together with its dead control surface:
  - no stage-7 `rotamer`, `placement_fixed_scalar`, or `placement_fixed_point_vector_only`;
  - no stage-7 `hybrid_sc_map` requirement;
  - no stage-7 SC-relaxation / proxy-control attributes that only served the old probabilistic proxy path.
- Revised decision after the 2026-04-02 cleanup pass:
  - remove the dead prep-only `hybrid_sc_map` export/validation path as well; the active workflow only needs `hybrid_control`, `hybrid_bb_map`, and `hybrid_env_topology`;
  - `src/box.cpp` and `src/box.h` were audited during this pass and left unchanged because they no longer contain the retired stage-7 proxy/rotamer subsystem, only active box/NPT infrastructure.
- `ALA` and `GLY` remain excluded from the new SC table because the trained library covers the current `18` non-empty canonical sidechain types; those residues still participate through the existing backbone/environment MARTINI path only.

### Execution Phases
- [x] Phase A: Build `parameters/ff_2.1/martini.h5` from the completed external `SC-training` results.
- [x] Phase B: Implement a new stage-7 SC table runtime node and wire its box-scaling behavior.
- [x] Phase C: Replace the production-side script injection in `run_sim_1rkl.sh` with stage-7 CB + SC-table setup and remove the legacy SC rotamer path.
- [x] Phase D: Update the production helper workflow to match the new stage-7 path.
- [x] Phase E: Remove the remaining parallel protein-MARTINI force path so stage 7 uses direct Upside-carried BB/SC interactions only.
- [x] Phase F: Verify `martini.h5` generation, build integrity, stage-7 HDF5 injection structure, and workflow/runtime behavior for the direct-Upside design.

### Known Errors / Blockers
- No code-path blocker remains for stage-7 SC-table injection in this environment:
  - the repaired repository `.venv` provides both `h5py` and `tables`;
  - the stage-7 injector now falls back to the protein ITP when fresh prepared files omit `/input/sequence`.
- Reduced-step full-workflow runs are still not a stability benchmark for the hybrid physics:
  - the intentionally shortened `6.4-6.6` relaxation settings used for fast verification can drive the pre-production MARTINI state to enormous energies before stage 7;
  - that does not block the stage-7 integration itself, but any physical assessment still needs the intended relaxation horizon / benchmark settings.

### Review
- 2026-04-02 legacy-code removal verification:
  - physically removed the dead probabilistic SC / rotamer / placement subsystem from `src/martini.cpp` instead of leaving it compiled behind runtime gates;
  - removed the dead prep-only `hybrid_sc_map` generation and validation path from `example/16.MARTINI/prepare_system_lib.py`, `example/16.MARTINI/prepare_system.py`, and `example/16.MARTINI/validate_hybrid_mapping.py`;
  - audited `src/box.cpp` and `src/box.h` and left them intact because the remaining code there is still-active barostat/box plumbing rather than disabled hybrid legacy logic;
  - verification:
    - `cmake --build obj` passed,
    - `python -m py_compile example/16.MARTINI/prepare_system_lib.py example/16.MARTINI/prepare_system.py example/16.MARTINI/validate_hybrid_mapping.py` passed under `.venv`,
    - `bash -n example/16.MARTINI/run_sim_1rkl.sh` and `bash -n example/16.MARTINI/test_prod_run_sim_1rkl.sh` passed,
    - a fresh prep run in `/tmp/hybrid_mapping_bb_only` produced a mapping file whose `/input` group contains only `hybrid_bb_map`, `hybrid_control`, and `hybrid_env_topology`,
    - a shortened active-workflow run in `/tmp/legacy_cleanup_short` accepted that mapping schema and completed the full current `run_sim_1rkl.sh` ladder through stage `7.0`.
- 2026-04-02 production-drift fix verification:
  - restored the startup protein-feedback ramp on the direct stage-7 force paths in `src/martini.cpp`:
    - direct `BB`/environment MARTINI forces now scale only the protein-side feedback by `compute_sc_backbone_feedback_mix(...)`,
    - `martini_sc_table_potential` now does the same for `CB` feedback while leaving the environment-side force unchanged;
  - removed the remaining active protein-internal MARTINI proxy-proxy path during stage 7 by making `allow_protein_pair_by_rule(...)` reject all protein-proxy pairs in the active hybrid runtime;
  - targeted replay verification from the same saved production handoff (`example/16.MARTINI/outputs/martini_test_1rkl_hybrid/checkpoints/1rkl.stage_7.0.up`) now stays negative through the user’s failure window:
    - archived current run at step `2050`: `martini_potential 26387.78`, `total 26216.96`;
    - patched replay at step `2050`: `martini_potential -25985.45`, `total -26147.28`.
- Implemented `example/16.MARTINI/build_sc_martini_h5.py` and generated `parameters/ff_2.1/martini.h5` from the completed external `SC-training` results.
- Implemented `example/16.MARTINI/inject_sc_table_stage7.py` so stage 7 now:
  - preserves the existing Upside backbone-node augmentation,
  - deletes the legacy SC production path,
  - injects deterministic CB coordinates and the new `martini_sc_table_potential`.
- Updated both `example/16.MARTINI/run_sim_1rkl.sh` and `example/16.MARTINI/test_prod_run_sim_1rkl.sh` to use the new stage-7-only SC table path and auto-build the library if only the assembled training JSON is available.
- Implemented the new `martini_sc_table_potential` runtime in `src/martini.cpp`:
  - CB-to-env radial table evaluation,
  - equal-and-opposite force accumulation,
  - additive backbone feedback through the CB coord node,
  - participation in NPT box scaling.
- Removed the remaining parallel protein proxy path from the active stage-7 runtime/workflow:
  - `martini_potential` now skips protein `SC` proxy MARTINI pairs entirely during active hybrid production;
  - protein `BB` proxy MARTINI pairs now use the Upside `CA` carrier coordinate as the interaction site and accumulate force directly onto `CA` plus the dry-MARTINI partner;
  - the active `hybrid_bb_map` is now CA-only in both exported mapping files and injected runtime stage files (`atom_mask=[0,1,0,0]`, `weights=[0,1,0,0]`), and runtime parsing now rejects any remaining active `N/C/O` carrier rows;
  - stage-file injection no longer copies `hybrid_sc_map` into runtime stages, the old shell-level SC control helper is gone from the active workflow path, and the legacy stage-7 rotamer augmenter now hard-fails if called.
- Verification completed in this environment:
  - `python3 -m py_compile` passed for the new Python scripts;
  - `bash -n` passed for both workflow shell scripts;
  - `cmake --build obj` passed;
  - `h5ls -r parameters/ff_2.1/martini.h5` confirmed the generated library shape.
- Regression revert verification completed after the user reported a broken minimization stage:
  - reverted the attempted native-mass switch and the CA-only active BB remap;
  - a fresh shortened `example/16.MARTINI/run_sim_1rkl.sh` run again reaches stage `6.0` with a normal minimization startup:
    - `Initial potential energy (Upside/MARTINI/Total): 0.00/140510.94/140510.94`
  - this replaces the broken regression signature reported by the user:
    - `Initial potential energy (Upside/MARTINI/Total): 0.00/30400945802214465601536.00/...`
- Native-mass-only restoration verification completed on 2026-04-02:
  - restored the native mass path while leaving the weighted 4-carrier BB map in place;
  - a fresh shortened `example/16.MARTINI/run_sim_1rkl.sh` run still starts stage `6.0` at:
    - `Initial potential energy (Upside/MARTINI/Total): 0.00/140510.94/140510.94`
  - fresh `1rkl.stage_6.0.up` in `/tmp/martini_mass_restore` again shows native dry-MARTINI mass values in `/input/mass` (`72.0` and `45.0` for the physical dry particles).
- Reduced-mass restoration verification completed after the user reported box blow-up with original masses:
  - restored reduced simulation masses (`ff_mass / 12`) in preparation and stage-file augmentation while leaving the scientific-notation logging change in place;
  - a fresh shortened `example/16.MARTINI/run_sim_1rkl.sh` run in `/tmp/martini_mass_reduced` again starts stage `6.0` at:
    - `Initial potential energy (Upside/MARTINI/Total): 0.00/140510.94/140510.94`
  - direct HDF5 inspection of `1rkl.stage_6.0.up` confirms the reduced mass values are back in `/input/mass`:
    - `4082` particles at `6.0`,
    - `8` particles at `3.75`.
- Logging/readability follow-up on 2026-04-02:
  - the huge `martini_potential` / `total` values seen in current stage-7 logs are real runtime values, not a `printf` type mismatch or summation bug;
  - `src/main.cpp` now formats large-magnitude energies in scientific notation for readability only.
- Additional end-to-end verification completed after the `.venv` repair:
  - the focused production helper `example/16.MARTINI/test_prod_run_sim_1rkl.sh` now runs end to end from a real `6.6 -> 7.0` handoff, injects the new SC table successfully, and completes a short stage-7 production replay;
  - direct HDF5 inspection confirmed stage gating:
    - `6.6` contains no `martini_sc_table_potential`, no `placement_fixed_point_only_CB`, no `hybrid_sc_map`, and keeps hybrid activation disabled;
    - `7.0.prepared` contains `martini_sc_table_potential` and `placement_fixed_point_only_CB`, remains free of legacy `rotamer` nodes, and activates the hybrid path only for production;
  - the real top-level workflow `example/16.MARTINI/run_sim_1rkl.sh` completes end to end under shortened verification settings with the direct-Upside stage-7 path.

## 2026-03-31 SC-training Workflow Addendum

### Project Goal
- Establish a new `SC-training/` workflow that derives dry-MARTINI sidechain-type to dry-MARTINI particle-type training data and assembled tables independently of the legacy SC back-mapping runtime path.
- Make the workflow runnable both locally and on Slurm.
- Use `example/16.MARTINI/run_sim_1rkl.sh` as the canonical benchmark workflow target once the trained tables are consumed by the hybrid runtime.

### Architecture & Key Decisions
- The new workflow will live under a dedicated top-level folder `SC-training/`.
- The workflow will treat the current SC back-mapping / rotamer-placement production path as legacy design context rather than a dependency.
- The replacement SC/dry-MARTINI interaction must be a two-way force term: each SC/dry pair contributes equal-and-opposite force to the dry particle and to the effective protein-side SC site, and the protein-side contribution must be projected back onto the Upside carriers and added to the existing backbone/sidechain forces before integration.
- Existing backbone/environment MARTINI coupling remains in scope and separate from the new SC table:
  - the effective protein backbone dry-MARTINI `BB`/`CA` anchor should still interact with surrounding non-protein dry-MARTINI particles;
  - the new SC table should not reintroduce generic protein SC-backbone MARTINI interactions, because protein internal SC-backbone coupling is already handled on the Upside side.
- Residue-sidechain dry-MARTINI bead definitions will be sourced from `example/16.MARTINI/martinize.py` using the `martini22` forcefield definitions already used by the example workflow.
- Dry-MARTINI pair parameters will be sourced from `example/16.MARTINI/ff_dry/dry_martini_v2.1*.itp`.
- The unit contract must be explicit and parameterized:
  - `SC-training/` artifacts remain in native dry-MARTINI units (`nm`, `kJ/mol`, `e`, dielectric `15`);
  - simulation-side conversion into Upside units must come from explicit parameters / attrs, not baked conversion numbers in training outputs or MARTINI runtime code;
  - the current MARTINI workflow expects `UPSIDE_MARTINI_ENERGY_CONVERSION` and `UPSIDE_MARTINI_LENGTH_CONVERSION` to be provided explicitly when generating simulation inputs.
- The initial training workflow will emphasize reproducible task generation, batch execution, and result assembly; any simplifying assumptions used for the first pass must be explicit in workflow outputs and README docs.
- Slurm support should follow the proven ConDiv workflow shape: explicit path resolution, environment setup, generated manifests/specs, array-friendly execution, and a submission wrapper rather than an ad hoc one-off `sbatch` command.

### Execution Phases
- [x] Phase A: Inspect and reuse the relevant MARTINI residue/type definitions and Slurm workflow patterns.
- [x] Phase B: Create `SC-training/` scripts for manifest generation, per-task training execution, and final result assembly.
- [x] Phase C: Add local and Slurm entrypoints that stage and run training jobs from the new folder.
- [x] Phase D: Add minimal documentation describing assumptions, outputs, and benchmark expectations.
- [x] Phase E: Verify script integrity and basic dry-run behavior.

### Known Errors / Blockers
- Runtime consumption of the new SC training table is not yet implemented in `src/martini.cpp`, so benchmark orchestration may need to stop at workflow staging/logging until that integration exists.
- That runtime integration still needs to enforce the new mechanical requirements explicitly:
  - SC/dry interactions must remain equal-and-opposite;
  - SC feedback must be added to existing Upside carrier gradients rather than replacing them;
  - backbone dry-`BB`/`CA` to environment dry-particle coupling must remain active while protein internal SC-backbone coupling stays excluded from the new table.
- The first-pass training workflow may need an explicit simplifying geometry assumption for residue-level sidechain reduction; that assumption must be documented, not hidden.
- `SC-training/` is now self-contained for training and Slurm staging, but the optional benchmark hook still targets `example/16.MARTINI/run_sim_1rkl.sh` and therefore still requires the full repository when benchmark execution is requested.

### Review
- Implemented a new `SC-training/` workflow with:
  - manifest generation from `example/16.MARTINI/martinize.py` (`martini22`) and `example/16.MARTINI/ff_dry/dry_martini_v2.1.itp`,
  - per-task residue/target training jobs,
  - assembled output tables,
  - local runner,
  - Slurm staging/submission wrappers,
  - benchmark staging against `example/16.MARTINI/run_sim_1rkl.sh`.
- Smoke verification succeeded in `/tmp/sc_training_smoke`:
  - manifest generation completed,
  - full local run completed,
  - assembled summary reported `684` tasks across `18` residues and `38` target particle types,
  - Slurm array, collector, and benchmark scripts were generated successfully with `--no-submit`.
- Follow-up workflow correction restored the default target set to the full dry-MARTINI atomtype list from the bundled forcefield after a temporary over-narrowing to non-ring types:
  - default table cardinality is `18 x 38 = 684` residue-target tasks;
  - manifest/task outputs now keep the full target set while retaining the explicit spherical surrounding-position samples around the sidechain center.
- Follow-up architecture review confirmed additional force-path requirements for the future runtime hookup:
  - the new SC/dry term must remain a two-way force and its protein-side feedback must be added on top of existing Upside forces;
  - backbone dry-`BB`/`CA` to environment dry-particle coupling remains required and should stay distinct from the new SC table;
  - dry-MARTINI unit handling is now defined as native-unit training plus parameterized simulation-side conversion, not a baked training/runtime numeric contract.
- Portability follow-up verification succeeded from a standalone folder copy under `/private/tmp`:
  - `SC-training/data/` now bundles the required dry-MARTINI nonbond parameter file and martini22 sidechain residue definitions;
  - `workflow.py init-run`, `workflow.py run-local`, `run_local.sh`, and stage-only `submit_remote_round.sh` all work from the copied folder without the parent repository;
  - generated Slurm array/collector scripts now reference the copied `SC-training/` folder directly instead of requiring `source.sh` or repo-root paths for training.
- Slurm-array verification succeeded from a fresh staged run:
  - `submit_remote_round.sh` now reaches `workflow.py submit-slurm` reliably even when repo-root `source.sh` would otherwise trip over an unset `MY_PYTHON`;
  - generated training script uses `#SBATCH --array=0-683` for the current manifest and dispatches one residue-target task per array index via `run-array-task --task-id "$SLURM_ARRAY_TASK_ID"`;
  - staged round manifest count and training manifest count matched exactly (`684` / `684`), with a separate collector script for post-array assembly.
- Batch-spooled launcher verification succeeded after the user-reported Midway failure:
  - `submit_remote_round.sh` and `run_local.sh` no longer trust `dirname "$0"` as the workflow root when executed from a Slurm spool copy;
  - the wrappers now locate the real `SC-training` directory by checking candidate directories for `workflow.py`, including `SLURM_SUBMIT_DIR` and `SLURM_SUBMIT_DIR/SC-training`;
  - simulated `sbatch submit_remote_round.sh` from inside `SC-training/` and simulated `sbatch SC-training/submit_remote_round.sh` from the repo root both staged cleanly and produced matching `684`-task manifests.
- Unit-contract follow-up verification succeeded:
  - training manifests now record native-unit policy only and no longer emit baked Upside conversion numbers;
  - `prepare_system_lib.py` now requires explicit conversion env vars and writes those MARTINI unit-conversion attrs into the simulation input;
  - `src/martini.cpp` now derives Coulomb scaling from those attrs instead of a literal hardcoded runtime constant.

### Phase 1: Geometric Mapping (Defining the Training Ground)
Before you can train the splines, you must define the physical relationship between the two models.
1. **Define the Upside Virtual Bead:** For each of the 20 standard amino acids, identify the Upside definition for the virtual side-chain bead position ($\vec{r}_{SC}$) and its normal vector ($\vec{n}$) relative to the backbone.
2. **Define the MARTINI Target:** For each amino acid, map the corresponding multi-bead Dry-MARTINI representation (e.g., mapping a Tryptophan side chain to its 4 constituent MARTINI beads). 

### Phase 2: Generating the 2D Energy Landscapes
You need to calculate the "ground truth" thermodynamic interactions using the Dry-MARTINI force field. Because Dry-MARTINI is an implicit solvent model, this can be done in a vacuum.
1. **Set Up the Grid:** Place your multi-bead MARTINI side chain at the origin. Place a single target Dry-MARTINI particle (e.g., a standard $C1$ lipid bead or $P4$ polar bead) in the space around it.
2. **Scan the Space:** Systematically move the single target bead across a 2D grid defined by:
   * **Distance ($r_{12}$):** The distance from the target bead to the Upside virtual center ($\vec{r}_{SC}$).
   * **Angle ($\theta_1$):** The angle between the Upside normal vector ($\vec{n}$) and the vector connecting the two beads ($\vec{n}_{12}$).
3. **Calculate Energy:** At every grid point, calculate the sum of the standard Dry-MARTINI Lennard-Jones and Coulomb interactions between the single target bead and the multi-bead side chain.
4. **Repeat:** Do this for all 20 amino acids against all relevant Dry-MARTINI particle types.

### Phase 3: Mathematical Spline Fitting
You must compress those 2D energy landscapes into the modified, factorized Upside equation.
1. **Set Up the Equation:** Use the simplified hybrid Hamiltonian that drops the second angle term:
   $$V_{hybrid} = V_{radial}(r_{12}) + ang_1(-n_1 \cdot n_{12}) \cdot V_{angular}(r_{12})$$
2. **Initialize Splines:** Set up $V_{radial}$, $ang_1$, and $V_{angular}$ as 1D natural cubic splines with adjustable control points (knots).
3. **Optimize:** Use a least-squares fitting algorithm or gradient descent to adjust the spline knots until the output of the $V_{hybrid}$ equation closely matches the 2D energy grid you generated in Phase 2.

### Phase 4: Engine Integration (C++ Implementation)
Append the new physics into the Upside MD engine's core loop.
1. **Introduce Real Particles:** Modify the Upside state to accept independent, mass-bearing Dry-MARTINI particles with standard $x, y, z$ coordinates and momenta.
2. **Calculate the Force:** In the force-evaluation loop, evaluate $V_{hybrid}$ and calculate its spatial gradient (force) and angular derivative (torque).
3. **Push the Real Bead:** Apply the negative spatial gradient directly to the Dry-MARTINI particle's momentum.
4. **Push the Ghost Bead:** Apply the exact inverse spatial vector to the Upside virtual side chain ($\vec{r}_{SC}$). Apply the angular derivative as a torque to the virtual normal vector ($\vec{n}$).
5. **Distribute to the Backbone:** Feed the virtual force and torque into Upside's existing chain-rule matrices to distribute the physical forces onto the $N$, $C\alpha$, and $C$ backbone atoms.

### Phase 5: Validation and Thermodynamic Tuning
Run test simulations to ensure the hybrid physics do not break the model.
1. **Check Momentum Conservation:** Monitor the total momentum of the simulation box. If the sum of forces on the MARTINI beads and the Upside backbone atoms does not equal exactly zero at every time step, there is a bug in your chain-rule implementation.
2. **Tune the Scaling Factor ($\kappa$):** Introduce a global scaling multiplier to the hybrid potential: $V_{total} = \kappa(V_{hybrid})$. Run a folded protein in a box of your Dry-MARTINI particles. If the protein denatures because the MARTINI interactions are overwhelmingly strong compared to Upside's internal folding potentials, lower $\kappa$ until the thermodynamic balance is restored.

What specific biological system or environment are you aiming to simulate first once this hybrid engine is built?
