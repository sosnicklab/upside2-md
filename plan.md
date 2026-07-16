# Project Goal

Fix the dryMARTINI hybrid force field so protein-facing CGL-target interactions
use the same physical single-lipid tail-compression state already used by
SC-CGL. Regenerate the installed DOPC force-field tables, rebuild fresh hybrid
stage-7 inputs, and validate the resulting 1RKL/1AFO stage-7.0 trajectories.

# Architecture & Key Decisions

- Keep the physical model unchanged except for restoring the missing
  protein-facing target-CGL tail-relief term. Do not change timestep,
  temperature, protein mass, damping, or add restraints.
- Scope the conservative fix to protein-facing `cg_lipid_target` only.
  `cg_lipid_target_base` remains the base table for non-protein particles.
- Use the existing source-conditioned target tail-relaxation machinery in
  `py/martini_build_tables.py`; do not introduce a new fitted membrane term or
  target-specific heuristic.
- Regenerated `dopc.h5` must preserve the current CGL-CGL and SC-CGL contracts,
  including the compressed single-CGL state and explicit hidden compaction
  coordinate.
- Stage injection must fail loudly if protein-facing `cg_lipid_target` would
  silently fall back to the base tensor while explicit CGL compaction is active.
- Validate the fix with fresh hybrid stage-7.0 runs for both 1RKL and 1AFO.

# Execution Phases

- [x] Update the base DOPC table builder so `cg_lipid_target` writes the
  single-CGL compaction payload and so the base SC builder preserves the
  compressed-state branch.
- [x] Tighten stage injection/validation so protein-facing `cg_lipid_target`
  requires the explicit compaction payload when `cgl_compaction_state` is in use.
- [x] Regenerate `parameters/dryMARTINI/dopc.h5`, backing up the current file,
  and verify the new `cg_lipid_target` metadata/datasets.
- [x] Regenerate fresh hybrid stage-7.0 inputs for 1RKL and 1AFO and verify
  that protein-facing `cg_lipid_target` takes arguments
  `compose_vector6d,pos,cgl_compaction_state`.
- [x] Run fresh stage-7.0 trajectories and measure tilt/depth drift against the
  July 16, 2026 rejected outputs.
- [x] Update `example/16.MARTINI/cg_lipid_potentials.tex` so it documents the
  q-dependent protein-facing CGL-target relief after the new tables pass.
- [x] Record the accepted results and remaining risks in `progress.md` and
  `findings.md`.

# Known Errors / Blockers

- The repaired `parameters/dryMARTINI/dopc.h5` now fixes the original target-CGL
  asymmetry, but Thursday, July 16, 2026 validation is mixed:
  1RKL is corrected cleanly, while 1AFO no longer drifts out of the bilayer yet
  still shows a chain-resolved tilt redistribution that is not fully resolved.
- Fresh stage-7.0 validation metrics:
  1RKL old `34.1 -> 7.8 deg`, `0.62 -> 10.44 A`;
  1RKL new `34.4 -> 31.4 deg`, `0.26 -> 1.03 A`.
  1AFO old whole-protein `6.0 -> 3.6 deg`, `0.02 -> 5.26 A`;
  1AFO new whole-protein `7.0 -> 25.6 deg`, `2.44 -> 0.06 A`.
- Chain-resolved 1AFO input helix tilts from the OPM-oriented AA PDB are
  `A=22.6 deg`, `B=8.7 deg`. The corrected stage-7.0 rerun ends near
  `A=12.2 deg`, `B=26.3 deg`, so the depth bug is fixed but one helix still
  rotates away from its initial orientation.

## Revised Decisions

- Treat the current `martini-dev` HEAD `ad8981b` and the installed July 15,
  2026 artifacts as the implementation baseline for this fix; earlier
  transport-only notes around `52f637e` are no longer sufficient.
- The missing target-CGL compaction payload is the first conservative defect to
  repair because it matches the observed ejection/verticalization mechanism and
  is already defined physically in the repo.
- Use an exact target-group rebuild of `cg_lipid_target` through the patched
  builder path to replace the installed `dopc.h5` atomically, rather than
  waiting for a full CG-CG table rebuild during force-field validation.

# Review

- Installed artifact `parameters/dryMARTINI/dopc.h5` now carries
  `cg_lipid_target/{delta_extended,delta_compact,delta_compressed}` and the
  single-CGL compaction metadata required by the runtime.
- Reinjected stage-7.0 prepared files confirm protein-facing
  `cg_lipid_target` now uses arguments
  `compose_vector6d,pos,cgl_compaction_state`, while non-protein targets remain
  on `cg_lipid_target_base`.
- Fresh validation outputs were written to:
  `example/16.MARTINI/outputs/martini_1rkl_hybrid_target_fix_20260716/1rkl.stage_7.0.vtf`
  and
  `example/16.MARTINI/outputs/martini_1afo_hybrid_target_fix_20260716/1afo.stage_7.0.vtf`.
