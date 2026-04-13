# Hybrid Interface Sweep Workflow

## Project Goal
- Create a root-level workflow that sweeps `PROTEIN_ENV_INTERFACE_SCALE` across the real hybrid `example/16.MARTINI/run_sim_1rkl.sh` path.
- Support local execution and Slurm array submission.
- Keep the implementation self-contained under `hybrid-interface-sweep/`.

## Architecture & Key Decisions
- Reuse `example/16.MARTINI/run_sim_1rkl.sh` directly instead of reimplementing the stage ladder.
- Run each `(interface_scale, replicate)` task in its own task-local `RUN_DIR`.
- Store a manifest plus per-task result JSONs and a compact assembled status summary.
- Follow the path-resolution and Slurm submission pattern from `bilayer-lateral-diffusion/`.
- Capture a whitelist of relevant `run_sim_1rkl.sh` env overrides at `init-run` time so shortened smoke settings and production settings are reproducible.

## Execution Phases
- [x] Phase 1: Create the workflow folder and local tracker files.
- [x] Phase 2: Implement `workflow.py` for manifest generation, task execution, aggregation, and Slurm staging.
- [x] Phase 3: Add `run_local.sh`, `submit_remote_round.sh`, and `README.md`.
- [x] Phase 4: Run static verification and a reduced smoke run.
- [x] Phase 5: Record review notes, results, and findings.

## Known Errors / Blockers
- No blocker is known yet, but the workflow must keep using the real hybrid shell entrypoint instead of drifting into a duplicated Python-only simulation path.

## Review
- Implemented:
  - `workflow.py`
  - `run_local.sh`
  - `submit_remote_round.sh`
  - `README.md`
  - local tracker files
- The workflow:
  - creates a manifest of `(interface_scale, replicate)` tasks,
  - runs `example/16.MARTINI/run_sim_1rkl.sh` directly,
  - writes one task-local `RUN_DIR` per task,
  - assembles compact execution summaries,
  - stages Slurm array and collector scripts.
- Verification:
  - `python3 -m py_compile hybrid-interface-sweep/workflow.py`
  - `bash -n hybrid-interface-sweep/run_local.sh`
  - `bash -n hybrid-interface-sweep/submit_remote_round.sh`
  - no-submit Slurm staging under `/tmp/hybrid_interface_sweep_smoke`
  - reduced local smoke run with one task at `interface_scale = 0.85`
  - confirmed the smoke task produced task-local `1rkl.stage_7.0.up` and a successful result JSON.
