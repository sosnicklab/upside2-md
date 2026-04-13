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

## Post-Run Analysis Extension

### Project Goal
- Add a post-run analysis path for completed sweep outputs.
- Measure protein and bilayer lateral diffusion signals from `stage_7.0.up`.
- Support local execution and Slurm submission for the analysis phase.

### Architecture & Key Decisions
- Discover completed `stage_7.0.up` files from the existing task tree.
- Analyze protein lateral COM MSD relative to bilayer COM using protein `CA` carrier atoms from `hybrid_bb_map`.
- Analyze lipid `PO4` lateral MSD as a bilayer fluidity guardrail.
- Keep the analysis post hoc and write its outputs under a separate analysis tree.

### Execution Phases
- [x] Phase 6: Extend `workflow.py` with analysis discovery, task execution, aggregation, and Slurm staging.
- [x] Phase 7: Add analysis wrapper/docs.
- [x] Phase 8: Run static verification and a reduced analysis smoke run.
- [x] Phase 9: Record review notes, results, and findings.

### Known Errors / Blockers
- The analysis must tolerate the current stage-file schema:
  - no `particle_class == LIPID` assumption,
  - no assumption that an output box dataset is present,
  - no assumption that protein carriers should be selected from proxy MARTINI atoms.

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
- Analysis extension implemented:
  - `workflow.py` now provides:
    - analysis manifest discovery,
    - per-task stage-7 MSD analysis,
    - assembled analysis CSV/JSON outputs,
    - analysis Slurm array and collector staging.
  - added:
    - `run_analysis_local.sh`
    - `submit_analysis.sh`
  - updated `README.md` with local and Slurm analysis usage.
- Analysis verification:
  - `python3 -m py_compile hybrid-interface-sweep/workflow.py`
  - `bash -n hybrid-interface-sweep/run_analysis_local.sh`
  - `bash -n hybrid-interface-sweep/submit_analysis.sh`
  - reduced real sweep run in `/tmp/hybrid_interface_sweep_analysis_smoke` with `PROD_70_NSTEPS=40`
  - `run_analysis_local.sh` created:
    - `analysis/analysis_manifest.json`
    - `analysis/results/tasks/scale0p85_r01.json`
    - `analysis/assembled/task_results.csv`
    - `analysis/assembled/condition_summary.csv`
    - `analysis/assembled/summary.json`
  - `submit_analysis.sh` staged:
    - `analysis/slurm/round_manifest.json`
    - `analysis/slurm/analyze_array.sbatch`
    - `analysis/slurm/collect_analysis.sbatch`
