#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import os
import shlex
import subprocess as sp
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Sequence


WORKFLOW_DIR = Path(os.path.abspath(__file__)).parent
REPO_ROOT = WORKFLOW_DIR.parent
HYBRID_RUN_SCRIPT = REPO_ROOT / "example" / "16.MARTINI" / "run_sim_1rkl.sh"

SCHEMA_MANIFEST = "hybrid_interface_sweep_manifest_v1"
SCHEMA_SLURM_ROUND = "hybrid_interface_sweep_slurm_round_v1"
SCHEMA_TASK_RESULT = "hybrid_interface_sweep_task_result_v1"

DEFAULT_INTERFACE_SCALES = [1.0, 0.85, 0.70, 0.55, 0.40, 0.25]
DEFAULT_REPLICATES = 3
DEFAULT_SEED = 20260413
DEFAULT_PDB_ID = "1rkl"

DEFAULT_PASSTHROUGH_ENV_KEYS = [
    "TEMPERATURE",
    "THERMOSTAT_TIMESCALE",
    "THERMOSTAT_INTERVAL",
    "MIN_60_MAX_ITER",
    "MIN_61_MAX_ITER",
    "EQ_62_NSTEPS",
    "EQ_63_NSTEPS",
    "EQ_64_NSTEPS",
    "EQ_65_NSTEPS",
    "EQ_66_NSTEPS",
    "PROD_70_NSTEPS",
    "EQ_TIME_STEP",
    "PROD_TIME_STEP",
    "MIN_TIME_STEP",
    "EQ_FRAME_STEPS",
    "PROD_FRAME_STEPS",
    "PROD_70_NPT_ENABLE",
    "PROD_70_BAROSTAT_TYPE",
    "PROD_70_BACKBONE_FIX_RIGID_ENABLE",
    "PRODUCTION_NONPROTEIN_HARD_SPHERE",
    "MARTINIZE_ENABLE",
    "MARTINIZE_FF",
    "PROTEIN_AA_PDB",
    "PROTEIN_CG_PDB",
    "PROTEIN_ITP",
    "BILAYER_PDB",
    "UNIVERSAL_PREP_MODE",
    "UPSIDE_MARTINI_ENERGY_CONVERSION",
    "UPSIDE_MARTINI_LENGTH_CONVERSION",
    "UPSIDE_EWALD_ENABLE",
    "UPSIDE_EWALD_ALPHA",
    "UPSIDE_EWALD_KMAX",
    "SC_ENV_LJ_FORCE_CAP",
    "SC_ENV_COUL_FORCE_CAP",
    "SC_ENV_RELAX_STEPS",
    "SC_ENV_BACKBONE_HOLD_STEPS",
    "SC_ENV_PO4_Z_HOLD_STEPS",
    "SC_ENV_PO4_Z_CLAMP_ENABLE",
    "SC_ENV_ENERGY_DUMP_ENABLE",
    "SC_ENV_ENERGY_DUMP_STRIDE",
    "SALT_MOLAR",
    "PROTEIN_LIPID_CUTOFF",
    "ION_CUTOFF",
    "BOX_PADDING_XY",
    "BOX_PADDING_Z",
    "PREP_SEED",
    "BB_AA_MIN_MATCHED_RESIDUES",
    "BB_AA_MAX_RIGID_RMSD",
    "SC_MARTINI_LIBRARY",
]


@dataclass
class Config:
    base_dir: Path
    pdb_id: str
    interface_scales: List[float]
    replicates: int
    seed: int
    passthrough_env: Dict[str, str]


def _now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def _manifest_path(base_dir: Path) -> Path:
    return base_dir / "sweep_manifest.json"


def _slurm_dir(base_dir: Path) -> Path:
    return base_dir / "slurm"


def _task_dir(base_dir: Path, code: str) -> Path:
    return base_dir / "tasks" / code


def _result_dir(base_dir: Path) -> Path:
    return base_dir / "results" / "tasks"


def _result_path(base_dir: Path, code: str) -> Path:
    return _result_dir(base_dir) / f"{code}.json"


def _assembled_dir(base_dir: Path) -> Path:
    return base_dir / "assembled"


def _write_json(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, sort_keys=True)
        fh.write("\n")


def _write_text(path: Path, content: str, executable: bool = False) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    if executable:
        path.chmod(0o755)


def _load_json(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def _load_manifest(base_dir: Path) -> Dict[str, Any]:
    path = _manifest_path(base_dir)
    if not path.exists():
        raise RuntimeError(f"Manifest not found: {path}. Run init-run first.")
    manifest = _load_json(path)
    if manifest.get("schema") != SCHEMA_MANIFEST:
        raise RuntimeError(f"Unexpected manifest schema in {path}: {manifest.get('schema')}")
    return manifest


def _float_list_arg(text: str) -> List[float]:
    values: List[float] = []
    for chunk in text.split(","):
        item = chunk.strip()
        if not item:
            continue
        values.append(float(item))
    if not values:
        raise argparse.ArgumentTypeError("Expected at least one numeric value")
    return values


def _format_float_tag(value: float, digits: int = 4) -> str:
    text = f"{value:.{digits}f}".rstrip("0").rstrip(".")
    if not text:
        text = "0"
    return text.replace("-", "m").replace(".", "p")


def _task_seed(base_seed: int, task_id: int) -> int:
    return int(base_seed + (task_id + 1) * 1009)


def _optional_sbatch_directives(phase: str) -> List[str]:
    directives: List[str] = []
    for key, flag in (
        ("PARTITION", "partition"),
        ("ACCOUNT", "account"),
        ("QOS", "qos"),
        ("CONSTRAINT", "constraint"),
        ("MEM", "mem"),
    ):
        value = os.environ.get(
            f"HYBRID_SWEEP_{phase}_{key}",
            os.environ.get(f"HYBRID_SWEEP_SBATCH_{key}", ""),
        ).strip()
        if value:
            directives.append(f"#SBATCH --{flag}={value}")
    return directives


def _shutil_which(binary: str) -> str | None:
    for path in os.environ.get("PATH", "").split(os.pathsep):
        if not path:
            continue
        candidate = Path(path) / binary
        if candidate.exists() and os.access(candidate, os.X_OK):
            return str(candidate)
    return None


def _run_sbatch(cmd: Sequence[str]) -> str:
    proc = sp.run(list(cmd), capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or proc.stdout.strip() or f"sbatch failed: {' '.join(cmd)}")
    return proc.stdout.strip().split(";", 1)[0]


def _prepend_env_path(env: Dict[str, str], key: str, value: str) -> None:
    existing = env.get(key, "")
    if existing:
        parts = [part for part in existing.split(os.pathsep) if part and part != value]
        env[key] = os.pathsep.join([value] + parts)
    else:
        env[key] = value


def _base_runtime_env() -> Dict[str, str]:
    env = os.environ.copy()
    env["UPSIDE_HOME"] = str(REPO_ROOT)
    env["PYTHONUNBUFFERED"] = "1"
    _prepend_env_path(env, "PATH", str(REPO_ROOT / "obj"))
    _prepend_env_path(env, "PYTHONPATH", str(REPO_ROOT / "py"))
    return env


def _call_logged(cmd: Sequence[str], log_path: Path, env: Dict[str, str]) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as fh:
        proc = sp.run(
            list(cmd),
            cwd=str(REPO_ROOT),
            env=env,
            stdout=fh,
            stderr=sp.STDOUT,
            text=True,
            check=False,
        )
    if proc.returncode != 0:
        detail = f"Command failed with exit code {proc.returncode}: {' '.join(cmd)}\nLog: {log_path}"
        try:
            tail_lines = log_path.read_text(encoding="utf-8", errors="replace").splitlines()[-40:]
        except OSError:
            tail_lines = []
        if tail_lines:
            detail = f"{detail}\nLast log lines:\n" + "\n".join(tail_lines)
        raise RuntimeError(detail)


def _normalize_env_key_list(values: Sequence[str]) -> List[str]:
    out: List[str] = []
    seen: set[str] = set()
    for raw in values:
        key = raw.strip()
        if not key or key in seen:
            continue
        seen.add(key)
        out.append(key)
    return out


def _capture_passthrough_env() -> Dict[str, str]:
    extra = _floatless_csv(os.environ.get("HYBRID_SWEEP_EXTRA_ENV_KEYS", ""))
    keys = _normalize_env_key_list(DEFAULT_PASSTHROUGH_ENV_KEYS + extra)
    out: Dict[str, str] = {}
    for key in keys:
        value = os.environ.get(key)
        if value is None or value == "":
            continue
        out[key] = value
    return out


def _floatless_csv(text: str) -> List[str]:
    items: List[str] = []
    for chunk in text.split(","):
        item = chunk.strip()
        if item:
            items.append(item)
    return items


def _build_config(args: argparse.Namespace, base_dir: Path) -> Config:
    scales = sorted(dict.fromkeys(float(x) for x in args.interface_scales), reverse=True)
    for scale in scales:
        if not (scale > 0.0):
            raise ValueError(f"interface scale must be > 0, got {scale}")
    replicates = int(args.replicates)
    if replicates < 1:
        raise ValueError(f"replicates must be >= 1, got {replicates}")
    return Config(
        base_dir=base_dir,
        pdb_id=str(args.pdb_id),
        interface_scales=scales,
        replicates=replicates,
        seed=int(args.seed),
        passthrough_env=_capture_passthrough_env(),
    )


def _build_manifest(config: Config) -> Dict[str, Any]:
    tasks: List[Dict[str, Any]] = []
    task_id = 0
    for interface_scale in config.interface_scales:
        for replicate in range(config.replicates):
            code = (
                f"scale{_format_float_tag(interface_scale)}_"
                f"r{replicate + 1:02d}"
            )
            tasks.append(
                {
                    "task_id": task_id,
                    "code": code,
                    "interface_scale": float(interface_scale),
                    "replicate": int(replicate + 1),
                    "seed": _task_seed(config.seed, task_id),
                }
            )
            task_id += 1
    return {
        "schema": SCHEMA_MANIFEST,
        "created_at_utc": _now_utc(),
        "base_dir": str(config.base_dir),
        "workflow_dir": str(WORKFLOW_DIR),
        "repo_root": str(REPO_ROOT),
        "settings": {
            "pdb_id": config.pdb_id,
            "interface_scales": config.interface_scales,
            "replicates": config.replicates,
            "seed": config.seed,
            "run_script": str(HYBRID_RUN_SCRIPT),
            "passthrough_env": config.passthrough_env,
        },
        "tasks": tasks,
    }


def _require_runtime_files() -> None:
    if not HYBRID_RUN_SCRIPT.exists():
        raise RuntimeError(f"Hybrid run script not found: {HYBRID_RUN_SCRIPT}")


def cmd_init_run(args: argparse.Namespace) -> int:
    _require_runtime_files()
    base_dir = Path(args.base_dir).expanduser().resolve()
    base_dir.mkdir(parents=True, exist_ok=True)
    config = _build_config(args, base_dir)
    manifest = _build_manifest(config)
    _write_json(_manifest_path(base_dir), manifest)
    print(f"Initialized hybrid interface sweep: {base_dir}")
    print(f"Task count: {len(manifest['tasks'])}")
    print(f"Run script: {manifest['settings']['run_script']}")
    return 0


def _run_task(base_dir: Path, manifest: Dict[str, Any], task: Dict[str, Any], overwrite: bool) -> Dict[str, Any]:
    result_path = _result_path(base_dir, task["code"])
    if result_path.exists() and not overwrite:
        existing = _load_json(result_path)
        if existing.get("success"):
            print(f"Skipping completed task {task['code']}")
            return existing

    settings = manifest["settings"]
    pdb_id = str(settings["pdb_id"])
    task_dir = _task_dir(base_dir, task["code"])
    run_dir = task_dir / "run"
    launcher_log = task_dir / "launcher.log"

    env = _base_runtime_env()
    env.update({str(k): str(v) for k, v in settings.get("passthrough_env", {}).items()})
    env["PDB_ID"] = pdb_id
    env["RUN_DIR"] = str(run_dir)
    env["CHECKPOINT_DIR"] = str(run_dir / "checkpoints")
    env["LOG_DIR"] = str(run_dir / "logs")
    env["HYBRID_PREP_DIR"] = str(run_dir / "hybrid_prep")
    env["PROTEIN_ENV_INTERFACE_SCALE"] = f"{float(task['interface_scale']):.8g}"
    env["SEED"] = str(int(task["seed"]))

    try:
        _call_logged(["bash", str(HYBRID_RUN_SCRIPT)], launcher_log, env)
        stage_70_file = run_dir / "checkpoints" / f"{pdb_id}.stage_7.0.up"
        if not stage_70_file.exists():
            raise RuntimeError(f"Missing expected stage-7 checkpoint: {stage_70_file}")
        result = {
            "schema": SCHEMA_TASK_RESULT,
            "created_at_utc": _now_utc(),
            "base_dir": str(base_dir),
            "task": task,
            "success": True,
            "run_dir": str(run_dir),
            "launcher_log": str(launcher_log),
            "stage_70_file": str(stage_70_file),
        }
        _write_json(result_path, result)
        return result
    except Exception as exc:
        result = {
            "schema": SCHEMA_TASK_RESULT,
            "created_at_utc": _now_utc(),
            "base_dir": str(base_dir),
            "task": task,
            "success": False,
            "run_dir": str(run_dir),
            "launcher_log": str(launcher_log),
            "error": str(exc),
        }
        _write_json(result_path, result)
        raise


def cmd_run_local(args: argparse.Namespace) -> int:
    base_dir = Path(args.base_dir).expanduser().resolve()
    manifest = _load_manifest(base_dir)
    tasks = manifest["tasks"]
    start_task = max(0, int(args.start_task))
    selected = tasks[start_task:]
    if args.max_tasks is not None:
        selected = selected[: max(0, int(args.max_tasks))]
    for task in selected:
        _run_task(base_dir, manifest, task, overwrite=bool(args.overwrite))
    if not args.no_assemble:
        assemble_results(base_dir)
    return 0


def cmd_run_array_task(args: argparse.Namespace) -> int:
    round_manifest = _load_json(Path(args.round_manifest).expanduser().resolve())
    if round_manifest.get("schema") != SCHEMA_SLURM_ROUND:
        raise RuntimeError(f"Unexpected round manifest schema: {round_manifest.get('schema')}")
    base_dir = Path(round_manifest["base_dir"]).expanduser().resolve()
    manifest = _load_manifest(base_dir)
    task_id = int(args.task_id)
    tasks = manifest["tasks"]
    if task_id < 0 or task_id >= len(tasks):
        raise RuntimeError(f"Task id out of range: {task_id} for {len(tasks)} tasks")
    _run_task(base_dir, manifest, tasks[task_id], overwrite=bool(args.overwrite))
    return 0


def _successful_results_from_dir(result_dir: Path) -> List[Dict[str, Any]]:
    results: List[Dict[str, Any]] = []
    if not result_dir.exists():
        return results
    for path in sorted(result_dir.glob("*.json")):
        payload = _load_json(path)
        if payload.get("success"):
            results.append(payload)
    return results


def assemble_results(base_dir: Path) -> int:
    manifest = _load_manifest(base_dir)
    results = _successful_results_from_dir(_result_dir(base_dir))
    assembled_dir = _assembled_dir(base_dir)
    assembled_dir.mkdir(parents=True, exist_ok=True)

    task_rows: List[Dict[str, Any]] = []
    grouped: Dict[float, List[Dict[str, Any]]] = {}
    for result in results:
        task = result["task"]
        row = {
            "task_id": int(task["task_id"]),
            "code": str(task["code"]),
            "interface_scale": float(task["interface_scale"]),
            "replicate": int(task["replicate"]),
            "seed": int(task["seed"]),
            "run_dir": str(result["run_dir"]),
            "launcher_log": str(result["launcher_log"]),
            "stage_70_file": str(result["stage_70_file"]),
        }
        task_rows.append(row)
        grouped.setdefault(row["interface_scale"], []).append(row)

    task_csv = assembled_dir / "task_results.csv"
    with task_csv.open("w", encoding="utf-8", newline="") as fh:
        fieldnames = [
            "task_id",
            "code",
            "interface_scale",
            "replicate",
            "seed",
            "run_dir",
            "launcher_log",
            "stage_70_file",
        ]
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in sorted(task_rows, key=lambda item: item["task_id"]):
            writer.writerow(row)

    summary_rows: List[Dict[str, Any]] = []
    expected_replicates = int(manifest["settings"]["replicates"])
    for scale, rows in sorted(grouped.items(), reverse=True):
        summary_rows.append(
            {
                "interface_scale": float(scale),
                "n_replicates_expected": expected_replicates,
                "n_replicates_completed": len(rows),
                "all_stage70_present": int(
                    all(Path(row["stage_70_file"]).exists() for row in rows)
                ),
            }
        )

    summary_csv = assembled_dir / "condition_summary.csv"
    with summary_csv.open("w", encoding="utf-8", newline="") as fh:
        fieldnames = [
            "interface_scale",
            "n_replicates_expected",
            "n_replicates_completed",
            "all_stage70_present",
        ]
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in summary_rows:
            writer.writerow(row)

    payload = {
        "created_at_utc": _now_utc(),
        "base_dir": str(base_dir),
        "n_tasks_total": int(len(manifest["tasks"])),
        "n_tasks_completed_successfully": len(task_rows),
        "n_scales_completed": len(summary_rows),
        "task_results_csv": str(task_csv),
        "condition_summary_csv": str(summary_csv),
    }
    _write_json(assembled_dir / "summary.json", payload)
    print(f"Assembled sweep results written under {assembled_dir}")
    return 0


def cmd_assemble_results(args: argparse.Namespace) -> int:
    base_dir = Path(args.base_dir).expanduser().resolve()
    return assemble_results(base_dir)


def _array_script_content(round_manifest_path: Path, base_dir: Path, n_tasks: int) -> str:
    output_path = _slurm_dir(base_dir) / "train-%A_%a.out"
    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name=hybrid-sweep-{base_dir.name}",
        f"#SBATCH --output={output_path}",
        f"#SBATCH --time={os.environ.get('HYBRID_SWEEP_TRAIN_WALLTIME', '24:00:00')}",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        f"#SBATCH --cpus-per-task={int(os.environ.get('HYBRID_SWEEP_CPUS_PER_TASK', '1'))}",
        f"#SBATCH --array=0-{n_tasks - 1}",
    ]
    lines.extend(_optional_sbatch_directives("TRAIN"))
    lines.extend(
        [
            "",
            "set -euo pipefail",
            "",
            f"SCRIPT_DIR={shlex.quote(str(WORKFLOW_DIR))}",
            f"ROUND_MANIFEST={shlex.quote(str(round_manifest_path))}",
            f"PROJECT_ROOT={shlex.quote(str(REPO_ROOT))}",
            'PYTHON_BIN=""',
            "",
            "if [ -f /etc/profile.d/modules.sh ]; then",
            "  source /etc/profile.d/modules.sh",
            "fi",
            "",
            "if command -v module >/dev/null 2>&1; then",
            "  module load python/3.11.9 || true",
            "  module load cmake || true",
            "  module load openmpi || true",
            '  module load "${HYBRID_SWEEP_HDF5_MODULE:-hdf5/1.14.3}" || true',
            "fi",
            "",
            'if [ -f "$SCRIPT_DIR/.venv/bin/activate" ]; then',
            '  source "$SCRIPT_DIR/.venv/bin/activate"',
            'elif [ -f "$PROJECT_ROOT/.venv/bin/activate" ]; then',
            '  source "$PROJECT_ROOT/.venv/bin/activate"',
            "fi",
            "",
            'if [ -n "${HYBRID_SWEEP_PYTHON:-}" ]; then',
            '  if command -v "${HYBRID_SWEEP_PYTHON}" >/dev/null 2>&1; then',
            '    PYTHON_BIN="$(command -v "${HYBRID_SWEEP_PYTHON}")"',
            "  else",
            '    PYTHON_BIN="${HYBRID_SWEEP_PYTHON}"',
            "  fi",
            'elif [ -n "${VIRTUAL_ENV:-}" ] && [ -x "${VIRTUAL_ENV}/bin/python3" ]; then',
            '  PYTHON_BIN="${VIRTUAL_ENV}/bin/python3"',
            "else",
            '  PYTHON_BIN="$(command -v python3)"',
            "fi",
            "",
            'if [ -z "${MY_PYTHON:-}" ]; then',
            '  py_path="$(command -v "$PYTHON_BIN" || true)"',
            '  if [ -n "$py_path" ]; then',
            '    export MY_PYTHON="$(cd "$(dirname "$py_path")/.." && pwd)"',
            "  fi",
            "fi",
            "set +u",
            'source "$PROJECT_ROOT/source.sh"',
            "set -u",
            'export PYTHONUNBUFFERED=1',
            'cmd=("$PYTHON_BIN" -u "$SCRIPT_DIR/workflow.py" run-array-task --round-manifest "$ROUND_MANIFEST" --task-id "$SLURM_ARRAY_TASK_ID")',
            'if [ "${HYBRID_SWEEP_OVERWRITE:-0}" = "1" ]; then',
            '  cmd+=(--overwrite)',
            "fi",
            '"${cmd[@]}"',
        ]
    )
    return "\n".join(lines) + "\n"


def _collector_script_content(base_dir: Path) -> str:
    output_path = _slurm_dir(base_dir) / "collect-%j.out"
    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name=hybrid-sweep-collect-{base_dir.name}",
        f"#SBATCH --output={output_path}",
        f"#SBATCH --time={os.environ.get('HYBRID_SWEEP_COLLECT_WALLTIME', '01:00:00')}",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        "#SBATCH --cpus-per-task=1",
    ]
    lines.extend(_optional_sbatch_directives("COLLECT"))
    lines.extend(
        [
            "",
            "set -euo pipefail",
            "",
            f"SCRIPT_DIR={shlex.quote(str(WORKFLOW_DIR))}",
            f"BASE_DIR={shlex.quote(str(base_dir))}",
            f"PROJECT_ROOT={shlex.quote(str(REPO_ROOT))}",
            'PYTHON_BIN=""',
            "",
            "if [ -f /etc/profile.d/modules.sh ]; then",
            "  source /etc/profile.d/modules.sh",
            "fi",
            "",
            "if command -v module >/dev/null 2>&1; then",
            "  module load python/3.11.9 || true",
            "  module load cmake || true",
            "  module load openmpi || true",
            '  module load "${HYBRID_SWEEP_HDF5_MODULE:-hdf5/1.14.3}" || true',
            "fi",
            "",
            'if [ -f "$SCRIPT_DIR/.venv/bin/activate" ]; then',
            '  source "$SCRIPT_DIR/.venv/bin/activate"',
            'elif [ -f "$PROJECT_ROOT/.venv/bin/activate" ]; then',
            '  source "$PROJECT_ROOT/.venv/bin/activate"',
            "fi",
            'if [ -n "${HYBRID_SWEEP_PYTHON:-}" ]; then',
            '  if command -v "${HYBRID_SWEEP_PYTHON}" >/dev/null 2>&1; then',
            '    PYTHON_BIN="$(command -v "${HYBRID_SWEEP_PYTHON}")"',
            "  else",
            '    PYTHON_BIN="${HYBRID_SWEEP_PYTHON}"',
            "  fi",
            'elif [ -n "${VIRTUAL_ENV:-}" ] && [ -x "${VIRTUAL_ENV}/bin/python3" ]; then',
            '  PYTHON_BIN="${VIRTUAL_ENV}/bin/python3"',
            "else",
            '  PYTHON_BIN="$(command -v python3)"',
            "fi",
            'export PYTHONUNBUFFERED=1',
            '"$PYTHON_BIN" -u "$SCRIPT_DIR/workflow.py" assemble-results --base-dir "$BASE_DIR"',
        ]
    )
    return "\n".join(lines) + "\n"


def cmd_submit_slurm(args: argparse.Namespace) -> int:
    base_dir = Path(args.base_dir).expanduser().resolve()
    manifest = _load_manifest(base_dir)
    round_manifest = {
        "schema": SCHEMA_SLURM_ROUND,
        "created_at_utc": _now_utc(),
        "base_dir": str(base_dir),
        "manifest_path": str(_manifest_path(base_dir)),
        "n_tasks": len(manifest["tasks"]),
    }
    round_manifest_path = _slurm_dir(base_dir) / "round_manifest.json"
    _write_json(round_manifest_path, round_manifest)

    array_script = _slurm_dir(base_dir) / "train_array.sbatch"
    collector_script = _slurm_dir(base_dir) / "collect_results.sbatch"
    _write_text(array_script, _array_script_content(round_manifest_path, base_dir, len(manifest["tasks"])), executable=True)
    _write_text(collector_script, _collector_script_content(base_dir), executable=True)

    print(f"Round manifest: {round_manifest_path}")
    print(f"Array script: {array_script}")
    print(f"Collector script: {collector_script}")

    if args.no_submit:
        print("Staged Slurm scripts only; no jobs submitted.")
        return 0

    if not _shutil_which("sbatch"):
        raise RuntimeError("sbatch not available in PATH")

    train_job = _run_sbatch(["sbatch", "--parsable", str(array_script)])
    collect_job = _run_sbatch(["sbatch", "--parsable", f"--dependency=afterok:{train_job}", str(collector_script)])
    _write_json(
        _slurm_dir(base_dir) / "submission.json",
        {
            "created_at_utc": _now_utc(),
            "train_job_id": train_job,
            "collect_job_id": collect_job,
        },
    )
    print(f"Submitted training array job: {train_job}")
    print(f"Submitted collector job: {collect_job}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Hybrid interface scale sweep workflow")
    sub = parser.add_subparsers(dest="command", required=True)

    p_init = sub.add_parser("init-run", help="Initialize a new hybrid interface scale sweep")
    p_init.add_argument("--base-dir", default=str(WORKFLOW_DIR / "runs" / "default"))
    p_init.add_argument("--pdb-id", default=DEFAULT_PDB_ID)
    p_init.add_argument("--interface-scales", type=_float_list_arg, default=DEFAULT_INTERFACE_SCALES)
    p_init.add_argument("--replicates", type=int, default=DEFAULT_REPLICATES)
    p_init.add_argument("--seed", type=int, default=DEFAULT_SEED)
    p_init.set_defaults(func=cmd_init_run)

    p_local = sub.add_parser("run-local", help="Run tasks locally")
    p_local.add_argument("--base-dir", default=str(WORKFLOW_DIR / "runs" / "default"))
    p_local.add_argument("--max-tasks", type=int, default=None)
    p_local.add_argument("--start-task", type=int, default=0)
    p_local.add_argument("--overwrite", action="store_true")
    p_local.add_argument("--no-assemble", action="store_true")
    p_local.set_defaults(func=cmd_run_local)

    p_array = sub.add_parser("run-array-task", help="Run a single Slurm array task")
    p_array.add_argument("--round-manifest", required=True)
    p_array.add_argument("--task-id", required=True, type=int)
    p_array.add_argument("--overwrite", action="store_true")
    p_array.set_defaults(func=cmd_run_array_task)

    p_assemble = sub.add_parser("assemble-results", help="Aggregate completed task results")
    p_assemble.add_argument("--base-dir", default=str(WORKFLOW_DIR / "runs" / "default"))
    p_assemble.set_defaults(func=cmd_assemble_results)

    p_submit = sub.add_parser("submit-slurm", help="Stage and optionally submit Slurm scripts")
    p_submit.add_argument("--base-dir", default=str(WORKFLOW_DIR / "runs" / "default"))
    p_submit.add_argument("--no-submit", action="store_true")
    p_submit.set_defaults(func=cmd_submit_slurm)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
