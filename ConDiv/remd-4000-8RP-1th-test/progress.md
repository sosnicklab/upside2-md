# Progress Log

- 2026-02-12: Created `task_plan.md`, `findings.md`, and `progress.md` per workflow requirements.
- 2026-02-12: Started failure investigation for worker crash on `1mg4` during local run.
- 2026-02-12: Investigated `run.output` and worker logs; identified mixed historical logs and confirmed one failure mode where non-venv `python3` lacks `numpy`.
- 2026-02-12: Compared `ConDiv.py` with `ConDiv_original.py`; found duplicated `__main__` mode dispatch in modernized script (behavioral drift).
- 2026-02-12: Updated `run_local.sh` to execute with explicit `venv/bin/python3` and print selected executable.
- 2026-02-12: Updated `ConDiv.py` to remove duplicated dispatch block and add worker log excerpt printing on missing divergence output.
- 2026-02-12: Synced `test_00/ConDiv.py` with `ConDiv.py` so checkpoint restart uses fixed worker script.
- 2026-02-12: Validation: `python3 -m py_compile ConDiv.py test_00/ConDiv.py` passed in venv; `bash -n run_local.sh` passed.
- 2026-02-12: Debugged new worker crash for `3wcq`; fixed `compute_divergence` namedtuple mutation bug (`contrast.sheet = ...`) by appending rows to `contrast.sheet` list, matching original mutable-list pattern.
- 2026-02-12: Synced updated `ConDiv.py` into `test_00/ConDiv.py` and re-ran syntax validation (`py_compile`) successfully.
- 2026-02-12: Added gradient L2 norm diagnostic in `backprop_deriv` debug print (`DEBUG: rot_grad norm=...`).
- 2026-02-12: `test_00/` is currently absent in workspace, so sync/compile was validated on top-level `ConDiv.py` only for this edit.
- 2026-02-12: Added `run_remote.sh` for Slurm usage with resume-from-latest-checkpoint semantics and configurable restart iterations (`$1` or `RUN_STEPS`, default 152).
- 2026-02-12: Made `run_remote.sh` executable and validated script syntax with `bash -n`.
- 2026-02-12: Added `README.md` documenting init flow, Slurm `run_remote.sh` submission/resume workflow, iteration overrides, and log locations.
- 2026-02-12: Updated `README.md` with required cluster module loads (`gcc/12.2.0`, `cmake`, `openmpi`) before Slurm submission.
