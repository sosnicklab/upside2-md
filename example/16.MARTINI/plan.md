# Hybrid Packing Calibration

## Project Goal
- Add a direct packing-analysis path for `example/16.MARTINI/` so membrane packing can be checked before further timescale tuning.
- Keep the production workflow behavior conservative:
  - do not silently flip production NPT defaults yet,
  - but make packing-oriented NPT tests easy to run and evaluate.

## Architecture & Key Decisions
- Analyze actual `stage_7.0.up` outputs instead of inferring packing from diffusion alone.
- Reuse the stage-file schema already produced by `run_sim_1rkl.sh`:
  - `output/pos`,
  - `output/time`,
  - `output/box` when present,
  - otherwise the static `martini_potential.{x_len,y_len,z_len}` attrs.
- Report lightweight membrane packing observables that can be computed directly from the current dry-MARTINI topology:
  - XY box area,
  - area per lipid from `PO4`-identified lipid molecules,
  - PO4-to-PO4 leaflet thickness,
  - tail-bond orientational order relative to the bilayer normal.
- Keep the new analysis as a standalone Python helper in `example/16.MARTINI/` rather than embedding it into the simulation workflow.
- Expose production NPT controls more clearly in the example wrapper/docs, but leave the shell default unchanged until packing comparisons justify a default flip.

## Execution Phases
- [x] Phase 1: Audit the current hybrid workflow and stage-7 output schema for packing-analysis inputs and production NPT controls.
- [x] Phase 2: Implement a standalone packing-analysis utility for `stage_7.0.up` outputs.
- [x] Phase 3: Update the example wrapper/docs so production NPT packing tests are easy to run intentionally.
- [ ] Phase 4: Compare current hybrid outputs against a dry-MARTINI packing reference and choose production NPT settings from measured APL/thickness/order metrics.

## Known Errors / Blockers
- No local dry-MARTINI packing reference artifact is stored inside this repository, so the new analyzer can quantify hybrid packing immediately but cannot, by itself, decide the final target values.
- `source .venv/bin/activate` may still be unreliable on this machine if the environment was moved; the repo-local `.venv/bin/python` remains usable and should be preferred for the new analyzer when needed.

## Review
- Implemented:
  - `example/16.MARTINI/analyze_packing.py`
  - explicit production-NPT wrapper knobs in `example/16.MARTINI/run.py`
  - packing-analysis usage and calibration guidance in `example/16.MARTINI/readme.md`
- Verification:
  - `.venv/bin/python -m py_compile example/16.MARTINI/analyze_packing.py example/16.MARTINI/run.py`
  - `bash -n example/16.MARTINI/run_sim_1rkl.sh`
  - `.venv/bin/python example/16.MARTINI/analyze_packing.py example/16.MARTINI/outputs/martini_test_1rkl_hybrid/checkpoints/1rkl.stage_7.0.up --json-out /tmp/packing_summary.json --csv-out /tmp/packing_timeseries.csv`
- Current baseline from the existing `martini_test_1rkl_hybrid` stage-7 output:
  - `APL = 89.01 A^2`
  - `PO4 thickness = 32.57 A`
  - `tail order = 0.238`
  - `output/box` is absent in that NVT output, so the analyzer correctly fell back to static `martini_potential` box attrs.
