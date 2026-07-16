# Progress

## Current Phase

- Patched `py/martini_build_tables.py` so the base DOPC builder writes the
  single-CGL compaction payload into `cg_lipid_target`, preserves the
  compressed-state SC payload, and records the runtime/reference compaction
  metadata needed by the stage runtime.
- Patched `py/martini_prepare_system_lib.py` so protein-facing
  `cg_lipid_target` now fails loudly if `cgl_compaction_state` is active but the
  installed `dopc.h5` lacks target `delta_*` datasets.
- Rebuilt the installed `parameters/dryMARTINI/dopc.h5` target group through
  the patched builder path and verified that `cg_lipid_target` now carries
  `delta_extended`, `delta_compact`, `delta_compressed`, and the explicit
  single-CGL compaction provenance.
- Re-injected copies of the existing 1RKL/1AFO `stage_7.0.prepared.up` files
  and verified that protein-facing `cg_lipid_target` now takes arguments
  `compose_vector6d,pos,cgl_compaction_state`, while
  `cg_lipid_target_base` remains `compose_vector6d,pos`.
- Ran fresh stage-7.0 validation from those reinjected prepared states into:
  `example/16.MARTINI/outputs/martini_1rkl_hybrid_target_fix_20260716/`
  and
  `example/16.MARTINI/outputs/martini_1afo_hybrid_target_fix_20260716/`.
- Validation outcome on Thursday, July 16, 2026:
  1RKL is corrected cleanly
  (`34.4 -> 31.4 deg`, `0.26 -> 1.03 A`);
  1AFO depth drift is corrected
  (`2.44 -> 0.06 A`, last-100 mean `0.41 A`)
  but the chain-resolved tilt remains mixed.

## Next

- Inspect the remaining 1AFO helix-tilt redistribution before treating the
  dryMARTINI interface as fully fixed.
