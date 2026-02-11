# Progress Log

- 2026-02-11: Initialized task files (`task_plan.md`, `findings.md`, `progress.md`).
- 2026-02-11: Reviewed bilayer output log and diffs; traced barostat sign error; fixed Berendsen scaling sign in `src/box.cpp`.
- 2026-02-11: Reviewed `example/16.MARTINI/bilayer_output_1.txt`; confirmed box dimensions change during NPT equilibration and stabilize during production.
- 2026-02-11: Created `example/16.MARTINI/run_sim_1ubq.sh` from the bilayer workflow and set default PDB ID to 1ubq.
- 2026-02-11: Ran `run_sim_1ubq.sh`; preparation failed because `pdb/1ubq_proa.itp` is missing, causing unknown protein atom mapping for MET BB.
- 2026-02-11: Added CRYST1 record to `example/16.MARTINI/pdb/1ubq.MARTINI.pdb` (copied from AA PDB).
- 2026-02-11: Re-ran `run_sim_1ubq.sh`; NPT equilibration hit NaN potential at step 20 and barostat expanded box.
- 2026-02-11: Updated `example/16.MARTINI/run_sim_1ubq.sh` to use smaller time step and longer thermostat timescale.
- 2026-02-11: Read CHARMM-GUI Gromacs .mdp files for 1ubq minimization/equilibration/production; noted minimization nsteps=3000, dt=0.020, Berendsen tau_p=5.0, compressibility=4.5e-5 bar^-1, and soft-core settings. Updated task_plan.md and findings.md accordingly.
- 2026-02-11: Updated `run_sim_1ubq.sh` to use Gromacs-aligned NPT settings (tau_p=5.0, isotropic, compressibility=2.1782) and increased minimization iterations to 3000. Re-ran workflow; NPT equilibration progressed past 100 steps without NaNs (log shows stable potentials/pressures).
- 2026-02-11: Re-ran `example/16.MARTINI/run_sim_1ubq.sh` to completion request; run is in progress (NPT equilibration step 80/2000 reached, no NaNs so far).
- 2026-02-11: Continued run in foreground; NPT equilibration progressing slowly (reached step 80/2000 with stable energies/pressures, no NaNs).

## 2026-02-11
- Added MARTINI backbone rigid hold registration and atom-name selection in `src/martini.cpp`.
- Added CLI switch `--martini-hold-backbone` and wiring in `src/main.cpp`.
- Stored `atom_names` in MARTINI input generation for backbone selection in `example/16.MARTINI/prepare_martini.py`.
- Enabled backbone hold for minimization and softened NPT equilibration stages in `example/16.MARTINI/run_sim_1ubq.sh`.
- Tests not run (not requested).
- Adjusted backbone selection to prefer explicit BB index dataset (`/input/martini_backbone_bb_indices`) sourced from PDB atom names in `example/16.MARTINI/prepare_martini.py`.
- Updated MARTINI rigid hold reader in `src/martini.cpp` to use BB indices when present, with atom-name fallback.
- Added H5 metadata arrays for atom roles, interaction type names, and per-atom particle class in `example/16.MARTINI/prepare_martini.py`.
- Updated MARTINI backbone selection to use /input/atom_roles as the fallback (instead of /input/atom_names), keeping BB indices as preferred source.
- Enforced MARTINI backbone hold requires /input/atom_roles or /input/martini_backbone_bb_indices when martini potential is present.
- Enforced presence of /input/type for MARTINI potential initialization.
- Removed /input/atom_roles generation; MARTINI backbone hold now relies solely on /input/martini_backbone_bb_indices.
- Restored /input/atom_roles from PDB atom roles and removed BB-index dataset usage; backbone hold now derives BB selection from atom roles.
- Switched atom_roles storage back to S4 dtype per request.
- Added per-atom molecule index dataset (/input/molecule_ids) with description and added residue_ids description in MARTINI input writer.
- Updated molecule grouping to keep protein chains intact (group by segid/chain), while non-proteins remain grouped by residue id; adjusted molecule validation and protein residue counting.
