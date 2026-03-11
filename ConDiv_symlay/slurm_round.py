#!/usr/bin/env python3
"""Stage and finalize separate-job Slurm rounds for ConDiv_symlay training."""

from __future__ import annotations

import argparse
import json
import os
import shlex
import shutil
import subprocess as sp
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import ConDiv_mem as cd


ROUND_SCHEMA = "condiv_symlay_slurm_round_v1"
WORKER_SCHEMA = "condiv_symlay_worker_spec_v1"
CONTINUE_EXIT_CODE = 10


def _now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def _phase_env(phase: str, key: str) -> str:
    return os.environ.get(f"CONDIV_{phase}_{key}", os.environ.get(f"CONDIV_SBATCH_{key}", "")).strip()


def _optional_sbatch_directives(phase: str) -> List[str]:
    directives: List[str] = []
    for key, flag in (
        ("PARTITION", "partition"),
        ("ACCOUNT", "account"),
        ("QOS", "qos"),
        ("CONSTRAINT", "constraint"),
        ("MEM", "mem"),
    ):
        value = _phase_env(phase, key)
        if value:
            directives.append(f"#SBATCH --{flag}={value}")
    return directives


def _optional_sbatch_args(phase: str) -> List[str]:
    args: List[str] = []
    for key, flag in (
        ("PARTITION", "partition"),
        ("ACCOUNT", "account"),
        ("QOS", "qos"),
        ("CONSTRAINT", "constraint"),
        ("MEM", "mem"),
    ):
        value = _phase_env(phase, key)
        if value:
            args.append(f"--{flag}={value}")
    return args


def _write_executable(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    path.chmod(0o755)


def _find_latest_checkpoint(base_dir: Path) -> Path:
    checkpoint = (base_dir / "initial_checkpoint.pkl").resolve()
    round_dirs = sorted(base_dir.glob("epoch_*_minibatch_*"), reverse=True)
    for round_dir in round_dirs:
        candidate = round_dir / "checkpoint.pkl"
        if candidate.exists():
            checkpoint = candidate.resolve()
            break
    if not checkpoint.exists():
        raise RuntimeError(f"Checkpoint not found under {base_dir}")
    return checkpoint


def _load_state(checkpoint: Path, base_dir: Path) -> dict:
    state = cd.load_checkpoint(str(checkpoint))
    state = cd._apply_runtime_config(state)
    state["base_dir"] = str(base_dir)
    layer_manifest_path = state.get("layer_manifest_path")
    if layer_manifest_path:
        os.environ.setdefault("CONDIV_SYMLAY_LAYER_MANIFEST", str(Path(layer_manifest_path).expanduser().resolve()))
    return state


def _round_dir_for_state(state: dict) -> Path:
    return Path(state["base_dir"]).expanduser().resolve() / (
        f"epoch_{state['epoch']:02d}_minibatch_{state['i_mb']:02d}"
    )


def _worker_spec_path(slurm_dir: Path, code: str) -> Path:
    return slurm_dir / "worker_specs" / f"{code}.json"


def _array_script_path(slurm_dir: Path) -> Path:
    return slurm_dir / "simulate_array.sbatch"


def _round_manifest_path(slurm_dir: Path) -> Path:
    return slurm_dir / "round_manifest.json"


def _submission_record_path(slurm_dir: Path) -> Path:
    return slurm_dir / "round_submission.json"


def _next_submission_path(round_dir: Path) -> Path:
    return round_dir / "next_round_submission.json"


def _build_worker_spec(state: dict, code: str, target, round_dir: Path, param_files: Dict[str, str]) -> dict:
    return {
        "schema": WORKER_SCHEMA,
        "generated_at_utc": _now_utc(),
        "code": code,
        "direc": str(round_dir),
        "fasta": target.fasta,
        "native_path": target.native_path,
        "chain_break": target.breakfile_path,
        "thickness": float(state["membrane_thickness"]),
        "nail_file": target.nail_file,
        "init_path": target.init_path,
        "n_res": int(target.n_res),
        "chi": target.chi,
        "param_files": dict(param_files),
        "sim_time": float(state["sim_time"]),
        "n_replica": int(state["n_replica"]),
        "omp_threads": int(state["omp_threads"]),
        "launch_mode": "local",
        "project_root": state["project_root"],
        "ff_dir": state["ff_dir"],
    }


def _array_script_content(
    script_dir: Path,
    project_root: Path,
    manifest_path: Path,
    round_dir: Path,
    n_task: int,
    omp_threads: int,
) -> str:
    output_path = round_dir / "slurm" / "simulate-%A_%a.out"
    directives = [
        "#!/bin/bash",
        f"#SBATCH --job-name=condiv-sim-{round_dir.name}",
        f"#SBATCH --output={output_path}",
        f"#SBATCH --time={os.environ.get('CONDIV_SIM_WALLTIME', '36:00:00')}",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        f"#SBATCH --cpus-per-task={int(omp_threads)}",
        f"#SBATCH --array=0-{int(n_task) - 1}",
    ]
    directives.extend(_optional_sbatch_directives("SIM"))
    directives.extend(
        [
            "",
            "set -euo pipefail",
            "",
            f"PROJECT_ROOT={shlex.quote(str(project_root))}",
            f"SCRIPT_DIR={shlex.quote(str(script_dir))}",
            "",
            "if [ -f /etc/profile.d/modules.sh ]; then",
            "  source /etc/profile.d/modules.sh",
            "fi",
            "",
            "if command -v module >/dev/null 2>&1; then",
            "  module load python/3.11.9",
            "  module load cmake",
            "  module load openmpi",
            "fi",
            "",
            'if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then',
            '  echo "ERROR: SLURM_ARRAY_TASK_ID is required." >&2',
            "  exit 1",
            "fi",
            "",
            f"export OMP_NUM_THREADS={int(omp_threads)}",
            "export OMP_PROC_BIND=${OMP_PROC_BIND:-close}",
            "export OMP_PLACES=${OMP_PLACES:-cores}",
            "export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}",
            "export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}",
            "export NUMEXPR_NUM_THREADS=${NUMEXPR_NUM_THREADS:-1}",
            "export VECLIB_MAXIMUM_THREADS=${VECLIB_MAXIMUM_THREADS:-1}",
            "export PYTHONUNBUFFERED=1",
            "",
            'VENV_ACTIVATE="$SCRIPT_DIR/venv/bin/activate"',
            'if [ ! -f "$VENV_ACTIVATE" ]; then',
            '  VENV_ACTIVATE="$PROJECT_ROOT/.venv/bin/activate"',
            "fi",
            'source "$VENV_ACTIVATE"',
            'export PYTHONPATH="${PYTHONPATH:-}"',
            'source "$PROJECT_ROOT/source.sh"',
            'export CONDIV_PROJECT_ROOT="$PROJECT_ROOT"',
            'export PYTHONPATH="$SCRIPT_DIR:$PROJECT_ROOT/py:$PROJECT_ROOT/obj:${PYTHONPATH:-}"',
            "",
            f'export CONDIV_ROUND_MANIFEST={shlex.quote(str(manifest_path))}',
            'python3 -u "$SCRIPT_DIR/slurm_round.py" run-array-task '
            '--round-manifest "$CONDIV_ROUND_MANIFEST" '
            '--task-id "$SLURM_ARRAY_TASK_ID"',
        ]
    )
    return "\n".join(directives) + "\n"


def _stage_round(
    checkpoint: Path,
    base_dir: Path,
    run_steps_remaining: int,
    resubmit_count: int,
) -> Tuple[Path, Path, dict]:
    state = _load_state(checkpoint, base_dir)
    round_dir = _round_dir_for_state(state)
    slurm_dir = round_dir / "slurm"
    slurm_dir.mkdir(parents=True, exist_ok=True)

    minibatch = state["minibatches"][state["i_mb"]]
    param_files = cd.prepare_minibatch_inputs(
        state,
        state["param"],
        state["init_param_files"],
        str(round_dir),
        state["solver"],
        state["sim_time"],
    )

    script_dir = Path(__file__).resolve().parent
    project_root = Path(state["project_root"]).expanduser().resolve()
    worker_specs: Dict[str, str] = {}
    task_order: List[str] = []
    for task_id, (code, target) in enumerate(minibatch):
        spec = _build_worker_spec(state, code, target, round_dir, param_files)
        spec["task_id"] = int(task_id)
        spec_path = _worker_spec_path(slurm_dir, code)
        spec_path.parent.mkdir(parents=True, exist_ok=True)
        spec_path.write_text(json.dumps(spec, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        worker_specs[code] = str(spec_path)
        task_order.append(code)

    manifest = {
        "schema": ROUND_SCHEMA,
        "generated_at_utc": _now_utc(),
        "checkpoint": str(checkpoint),
        "base_dir": str(base_dir),
        "round_dir": str(round_dir),
        "epoch": int(state["epoch"]),
        "minibatch": int(state["i_mb"]),
        "run_steps_remaining": int(run_steps_remaining),
        "resubmit_count": int(resubmit_count),
        "script_dir": str(script_dir),
        "project_root": str(project_root),
        "layer_manifest_path": state.get("layer_manifest_path", str(base_dir / "layer_manifest.json")),
        "progress_log": os.environ.get("PROGRESS_LOG_JSONL", str(base_dir / "training_progress.jsonl")),
        "status_json": os.environ.get("STATUS_JSON", str(base_dir / "training_status.json")),
        "grad_threshold": float(os.environ.get("CONDIV_CONVERGENCE_GRAD_NORM", "1.0")),
        "update_threshold": float(os.environ.get("CONDIV_CONVERGENCE_UPDATE_NORM", "0.05")),
        "patience": int(os.environ.get("CONDIV_CONVERGENCE_PATIENCE", "3")),
        "auto_resubmit": 1 if os.environ.get("CONDIV_AUTO_RESUBMIT", "1") == "1" else 0,
        "auto_max_submissions": int(os.environ.get("CONDIV_AUTO_MAX_SUBMISSIONS", "0")),
        "update_walltime": os.environ.get("CONDIV_UPDATE_WALLTIME", "02:00:00"),
        "worker_specs": worker_specs,
        "task_order": task_order,
        "n_array_tasks": int(len(task_order)),
        "target_names": [code for code, _target in minibatch],
    }

    manifest_path = _round_manifest_path(slurm_dir)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    array_script_path = _array_script_path(slurm_dir)
    _write_executable(
        array_script_path,
        _array_script_content(
            script_dir=script_dir,
            project_root=project_root,
            manifest_path=manifest_path,
            round_dir=round_dir,
            n_task=len(task_order),
            omp_threads=int(state["omp_threads"]),
        ),
    )
    return manifest_path, array_script_path, manifest


def _run_sbatch(cmd: List[str]) -> str:
    proc = sp.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or proc.stdout.strip() or f"sbatch failed: {' '.join(cmd)}")
    return proc.stdout.strip().split(";", 1)[0]


def submit_round(
    base_dir: Path,
    checkpoint: Path | None,
    run_steps_remaining: int,
    resubmit_count: int,
    no_submit: bool,
) -> int:
    base_dir = base_dir.expanduser().resolve()
    checkpoint = checkpoint.expanduser().resolve() if checkpoint else _find_latest_checkpoint(base_dir)
    manifest_path, array_script_path, manifest = _stage_round(checkpoint, base_dir, run_steps_remaining, resubmit_count)

    print("Round manifest:", manifest_path)
    print("Round dir:", manifest["round_dir"])
    print("Checkpoint:", checkpoint)
    print("Targets:", ", ".join(manifest["target_names"]))
    print("Simulation array script:", array_script_path)

    if no_submit:
        print("Staged Slurm scripts only; no jobs submitted.")
        return 0

    if int(manifest["n_array_tasks"]) <= 0:
        raise RuntimeError("No array tasks were generated")
    if not shutil.which("sbatch"):
        raise RuntimeError("sbatch is not available in PATH")

    array_job_id = _run_sbatch(["sbatch", "--parsable", str(array_script_path)])
    print(f"Submitted simulation array job: {array_job_id}")

    dependency = f"afterany:{array_job_id}"
    update_output = str(Path(manifest["round_dir"]) / "slurm" / "update-%j.out")
    update_script = Path(manifest["script_dir"]) / "run_remote_update.sh"
    update_cmd = [
        "sbatch",
        "--parsable",
        f"--dependency={dependency}",
        f"--job-name=condiv-upd-e{manifest['epoch']:02d}-m{manifest['minibatch']:02d}",
        f"--output={update_output}",
        f"--time={manifest['update_walltime']}",
        "--nodes=1",
        "--ntasks=1",
        "--cpus-per-task=1",
        *_optional_sbatch_args("UPDATE"),
        (
            "--export=ALL,"
            f"CONDIV_PROJECT_ROOT={manifest['project_root']},"
            f"CONDIV_SYMLAY_LAYER_MANIFEST={manifest['layer_manifest_path']},"
            f"CONDIV_ROUND_MANIFEST={manifest_path}"
        ),
        str(update_script),
    ]
    update_job_id = _run_sbatch(update_cmd)
    print(f"Submitted update job: {update_job_id}")

    record = {
        "schema": ROUND_SCHEMA,
        "submitted_at_utc": _now_utc(),
        "checkpoint": str(checkpoint),
        "round_manifest": str(manifest_path),
        "simulation_array_job": {
            "job_id": array_job_id,
            "script": str(array_script_path),
            "n_array_tasks": int(manifest["n_array_tasks"]),
        },
        "update_job": {"job_id": update_job_id, "script": str(update_script), "dependency": dependency},
    }
    _submission_record_path(Path(manifest["round_dir"]) / "slurm").write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return 0


def _worker_argv_from_spec(spec: dict) -> List[str]:
    return [
        str(Path(cd.__file__).resolve()),
        "worker",
        spec["code"],
        spec["direc"],
        spec["fasta"],
        spec["native_path"],
        spec["chain_break"],
        str(spec["thickness"]),
        spec["nail_file"],
        spec["init_path"],
        str(spec["n_res"]),
        spec["chi"],
        cd.cp.dumps(spec["param_files"], protocol=4).hex(),
        str(spec["sim_time"]),
        str(spec["n_replica"]),
        str(spec["omp_threads"]),
        str(spec["launch_mode"]),
        spec["project_root"],
        spec["ff_dir"],
    ]


def run_worker(worker_spec: Path) -> int:
    spec = json.loads(worker_spec.read_text(encoding="utf-8"))
    if spec.get("schema") != WORKER_SCHEMA:
        raise RuntimeError(f"Unexpected worker spec schema in {worker_spec}")

    argv = _worker_argv_from_spec(spec)
    old_argv = sys.argv[:]
    try:
        sys.argv = argv
        cd.main_worker()
    finally:
        sys.argv = old_argv
    return 0


def run_array_task(round_manifest: Path, task_id: int) -> int:
    manifest = json.loads(round_manifest.read_text(encoding="utf-8"))
    if manifest.get("schema") != ROUND_SCHEMA:
        raise RuntimeError(f"Unexpected round manifest schema in {round_manifest}")

    task_order = list(manifest.get("task_order", []))
    if task_id < 0 or task_id >= len(task_order):
        raise RuntimeError(f"Task id {task_id} is out of range for {round_manifest}")

    code = task_order[task_id]
    worker_specs = manifest.get("worker_specs", {})
    worker_spec = worker_specs.get(code)
    if not worker_spec:
        raise RuntimeError(f"No worker spec found for task {task_id} ({code}) in {round_manifest}")

    print(f"Running array task {task_id} for {code}")
    return run_worker(Path(worker_spec).expanduser().resolve())


def _run_training_control(manifest: dict, checkpoint: Path) -> int:
    training_control = Path(manifest["script_dir"]) / "training_control.py"
    cmd = [
        sys.executable,
        str(training_control),
        "--base-dir",
        manifest["base_dir"],
        "--checkpoint",
        str(checkpoint),
        "--progress-log",
        manifest["progress_log"],
        "--status-json",
        manifest["status_json"],
        "--grad-threshold",
        str(manifest["grad_threshold"]),
        "--update-threshold",
        str(manifest["update_threshold"]),
        "--patience",
        str(manifest["patience"]),
        "--run-steps",
        "1",
        "--slurm-job-id",
        os.environ.get("SLURM_JOB_ID", ""),
        "--resubmit-count",
        str(manifest["resubmit_count"]),
    ]
    proc = sp.run(cmd, check=False)
    return int(proc.returncode)


def finalize_round(manifest_path: Path) -> int:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != ROUND_SCHEMA:
        raise RuntimeError(f"Unexpected round manifest schema in {manifest_path}")

    checkpoint = Path(manifest["checkpoint"]).expanduser().resolve()
    base_dir = Path(manifest["base_dir"]).expanduser().resolve()
    round_dir = Path(manifest["round_dir"]).expanduser().resolve()
    state = _load_state(checkpoint, base_dir)
    expected_round_dir = _round_dir_for_state(state)
    if expected_round_dir != round_dir:
        raise RuntimeError(
            f"Round manifest {round_dir} does not match checkpoint state {expected_round_dir}"
        )

    minibatch = state["minibatches"][state["i_mb"]]
    new_param, grad_stats = cd.finalize_minibatch_from_outputs(
        state,
        state["param"],
        state["init_param_files"],
        str(round_dir),
        minibatch,
        state["solver"],
    )
    state = cd.advance_training_state(state, str(round_dir), new_param, grad_stats)
    new_checkpoint = round_dir / "checkpoint.pkl"
    cd.save_checkpoint(str(new_checkpoint), state)
    print("Saved checkpoint:", new_checkpoint)

    control_rc = _run_training_control(manifest, new_checkpoint)
    if control_rc == 0:
        print("Training converged; no new round submitted.")
        return 0
    if control_rc != CONTINUE_EXIT_CODE:
        print(f"Training control failed with exit code {control_rc}")
        return control_rc

    next_submission = _next_submission_path(round_dir)
    if next_submission.exists():
        print("Next round already submitted:", next_submission)
        return 0

    if int(manifest["run_steps_remaining"]) == 1:
        print("Reached run-steps limit for this submission chain; stopping before convergence.")
        return 0

    if int(manifest["auto_resubmit"]) != 1:
        print("Auto-resubmission disabled; stopping before convergence.")
        return 0

    next_resubmit_count = int(manifest["resubmit_count"]) + 1
    auto_max_submissions = int(manifest["auto_max_submissions"])
    if auto_max_submissions > 0 and next_resubmit_count > auto_max_submissions:
        print(f"Reached CONDIV_AUTO_MAX_SUBMISSIONS={auto_max_submissions}; stopping before convergence.")
        return 0

    next_run_steps = int(manifest["run_steps_remaining"])
    if next_run_steps > 0:
        next_run_steps -= 1

    rc = submit_round(
        base_dir=base_dir,
        checkpoint=new_checkpoint,
        run_steps_remaining=next_run_steps,
        resubmit_count=next_resubmit_count,
        no_submit=False,
    )
    if rc == 0:
        next_submission.write_text(
            json.dumps(
                {
                    "submitted_at_utc": _now_utc(),
                    "from_round_manifest": str(manifest_path),
                    "from_checkpoint": str(new_checkpoint),
                    "next_resubmit_count": next_resubmit_count,
                    "next_run_steps_remaining": next_run_steps,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
    return rc


def main() -> int:
    parser = argparse.ArgumentParser(description="ConDiv_symlay separate-job Slurm round orchestration")
    subparsers = parser.add_subparsers(dest="command", required=True)

    submit_parser = subparsers.add_parser("submit-round", help="Stage one minibatch round and submit its jobs")
    submit_parser.add_argument("--base-dir", required=True, help="ConDiv run directory")
    submit_parser.add_argument("--checkpoint", help="Checkpoint to resume from; defaults to latest")
    submit_parser.add_argument(
        "--run-steps",
        type=int,
        default=0,
        help="Remaining rounds to auto-submit in this chain; 0 means no chain limit",
    )
    submit_parser.add_argument("--resubmit-count", type=int, default=0, help="Current auto-resubmission count")
    submit_parser.add_argument("--no-submit", action="store_true", help="Only generate scripts and manifest")

    worker_parser = subparsers.add_parser("run-worker", help="Run one simulation from a generated worker spec")
    worker_parser.add_argument("--worker-spec", required=True, help="Path to a generated worker spec JSON")

    array_task_parser = subparsers.add_parser("run-array-task", help="Run one array task from a round manifest")
    array_task_parser.add_argument("--round-manifest", required=True, help="Path to a generated round manifest JSON")
    array_task_parser.add_argument("--task-id", required=True, type=int, help="Array task id")

    finalize_parser = subparsers.add_parser("finalize-round", help="Finalize one completed round")
    finalize_parser.add_argument("--round-manifest", required=True, help="Path to a generated round manifest JSON")

    args = parser.parse_args()
    if args.command == "submit-round":
        checkpoint = Path(args.checkpoint) if args.checkpoint else None
        return submit_round(
            base_dir=Path(args.base_dir),
            checkpoint=checkpoint,
            run_steps_remaining=int(args.run_steps),
            resubmit_count=int(args.resubmit_count),
            no_submit=bool(args.no_submit),
        )
    if args.command == "run-worker":
        return run_worker(Path(args.worker_spec).expanduser().resolve())
    if args.command == "run-array-task":
        return run_array_task(
            Path(args.round_manifest).expanduser().resolve(),
            int(args.task_id),
        )
    if args.command == "finalize-round":
        return finalize_round(Path(args.round_manifest).expanduser().resolve())
    raise RuntimeError(f"Unsupported command: {args.command}")


if __name__ == "__main__":
    sys.exit(main())
