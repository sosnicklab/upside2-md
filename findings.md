# Findings

## Root Cause

- On Thursday, July 16, 2026, the current rejected production outputs are:
  `example/16.MARTINI/outputs/martini_1rkl_hybrid/1rkl.stage_7.0.vtf` and
  `example/16.MARTINI/outputs/martini_1afo_hybrid/1afo.stage_7.0.vtf`.
- Direct measurement from those trajectories shows:
  1RKL tilt `34.1 -> 6.9 deg`, protein depth `0.9 -> 10.7 A`;
  1AFO chain tilts remain shallower but protein depth still drifts upward.
- The current `.up` stage files wire `cg_lipid_rotamer_sc` with
  `cgl_compaction_state`, but protein-facing `cg_lipid_target` is injected with
  only `compose_vector6d,pos` and no target `delta_*` datasets.
- Installed `parameters/dryMARTINI/dopc.h5` therefore has a conservative
  asymmetry:
  `cg_lipid_sc` contains `delta_extended`, `delta_compact`,
  `delta_compressed`, and single-CGL compaction provenance, while
  `cg_lipid_target` contains only the base interaction tensor.
- The repo already contains a physical target tail-relaxation implementation:
  `_target_conditioned_tail_relaxation_effective_energy`,
  `_build_target_tail_relaxation_compaction_payload`, explicit target-q runtime
  support in `src/martini_cg_lipid.cpp`, and injection logic that will attach
  `cgl_compaction_state` when those datasets exist.
- `build_dopc_h5()` currently writes the base target tensor directly and does
  not apply the later single-CGL compaction retrofit path. The base builder is
  therefore the artifact-generation defect that must be fixed.
- Git history between `52f637e` and current HEAD `ad8981b` shows no recent
  target-compaction regression in the C++ runtime; the only relevant post-
  `52f637e` code change in this area was a July 14, 2026 preparation update for
  protein mass and CGL rotational thermostat defaults.

## Physical Interpretation

- The current bug is consistent with protein-facing lipids being allowed to
  relax away from SC overlap but not from backbone/target overlap.
- Restoring q-dependent target relief is a physical conservative-force repair,
  not a transport or thermostat change.
- The fix should apply only to protein-facing target particles; non-protein
  targets should remain on `cg_lipid_target_base`.

## Repair Result

- The installed `parameters/dryMARTINI/dopc.h5` now contains
  `cg_lipid_target/{delta_extended,delta_compact,delta_compressed}` plus the
  explicit single-CGL compaction attributes used by the runtime.  The repaired
  target group was generated through the patched builder path and swapped into
  the installed H5 atomically.
- Re-injected stage-7.0 prepared files for both 1RKL and 1AFO now wire
  protein-facing `cg_lipid_target` with
  `compose_vector6d,pos,cgl_compaction_state`. Non-protein targets stay on the
  base tensor through `cg_lipid_target_base`.
- Fresh Thursday, July 16, 2026 validation from the reinjected prepared states
  shows that the original 1RKL failure mode is corrected:
  old `34.1 -> 7.8 deg`, `0.62 -> 10.44 A`;
  new `34.4 -> 31.4 deg`, `0.26 -> 1.03 A`.
- Fresh Thursday, July 16, 2026 validation shows that 1AFO no longer walks out
  of the bilayer:
  old whole-protein depth `0.02 -> 5.26 A`;
  new whole-protein depth `2.44 -> 0.06 A` with last-100 mean `0.41 A`.
- The 1AFO tilt outcome is still mixed.  Whole-protein PCA changes from
  `7.0 -> 25.6 deg` in the corrected run, but that metric hides the dimer
  asymmetry.  Chain-resolved helix tilts from the input AA structure are
  `A=22.6 deg`, `B=8.7 deg`; the corrected final state ends near
  `A=12.2 deg`, `B=26.3 deg`.  The depth bug is fixed, but a chain-tilt issue
  remains.

## Lessons

- When a hidden-state conservative term exists in both builder and runtime code,
  verify that the installed H5 artifact actually carries the required datasets
  before trusting trajectory behavior.
- Do not accept a stage injection path that silently downgrades a protein-facing
  interaction from explicit hidden-state coupling to a base table when the
  physical compaction state is active.
- For dimeric membrane proteins such as 1AFO, whole-protein PCA can hide
  chain-specific tilt redistribution. Validate both bilayer depth and
  chain-resolved helix tilt before calling an interface fix complete.
