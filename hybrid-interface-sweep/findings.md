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

## Lessons
- When a new sweep is supposed to exercise an existing scientific workflow, orchestrate the real entrypoint rather than re-encoding the stage logic in a second place.
