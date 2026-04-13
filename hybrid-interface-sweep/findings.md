# Findings

## 2026-04-13 (Initial Design)
- The correct execution surface for the sweep is `example/16.MARTINI/run_sim_1rkl.sh`, because that script already owns:
  - hybrid preparation,
  - stage-file generation,
  - stage `7.0` activation,
  - `PROTEIN_ENV_INTERFACE_SCALE`,
  - trajectory extraction.
- The right Slurm pattern already exists in `bilayer-lateral-diffusion/`:
  - a manifest created at `init-run`,
  - one array task per manifest row,
  - a collector job that assembles per-task results after the array finishes.

## 2026-04-13 (Implementation Details)
- The new workflow can stay thin because `run_sim_1rkl.sh` is already fully env-configurable for:
  - `PROTEIN_ENV_INTERFACE_SCALE`
  - `RUN_DIR`
  - simulation-length overrides
  - thermostat and production settings.
- Capturing a whitelist of relevant env overrides at `init-run` time is enough to make shortened smoke runs and production reruns reproducible without turning the new workflow into a second copy of the hybrid parameter surface.
- The reduced smoke run verified that the workflow can drive the canonical hybrid shell path all the way to a real task-local `stage_7.0.up` checkpoint.

## 2026-04-13 (Analysis Design)
- The active hybrid stage-7 outputs need schema-aware analysis:
  - protein carrier atoms are best selected from `hybrid_bb_map` rather than proxy MARTINI roles,
  - lipid membership is recoverable from non-protein `PO4` beads plus `molecule_ids`,
  - box lengths may need to fall back to `input/potential/martini_potential` attrs when `/output/box` is absent.
- A practical first analysis target for the sweep is diffusion, not a direct viscosity estimator:
  - protein lateral COM diffusion relative to the bilayer COM is the main sweep signal,
  - lipid `PO4` lateral diffusion is the bilayer guardrail,
  - the analysis can be assembled into per-task and per-condition CSV summaries for later offline comparison.
- The external `/Users/yinhan/Documents/bilayer-lateral-diffusion` workflow uses the intended acceleration pattern:
  - a manifest/init step,
  - one Slurm array element per task,
  - a dependent collector job after the array.
- The hybrid sweep now mirrors that pattern explicitly for both simulation and analysis submission.

## Lessons
- When a new sweep is supposed to exercise an existing scientific workflow, orchestrate the real entrypoint rather than re-encoding the stage logic in a second place.
- When adding post-run analysis, derive atom selection and box handling from the actual generated stage files before writing the analyzer.
- When the user points to a reference Slurm workflow, inspect its live array-submission path and match that pattern explicitly rather than relying on informal similarity.
