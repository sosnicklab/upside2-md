#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shlex
import shutil
import subprocess as sp
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence

import h5py
import numpy as np


WORKFLOW_DIR = Path(os.path.abspath(__file__)).parent
REPO_ROOT = WORKFLOW_DIR.parent
PREP_SCRIPT = REPO_ROOT / "py" / "martini_prepare_system.py"
UPSIDE_EXECUTABLE = REPO_ROOT / "obj" / "upside"

SCHEMA_MANIFEST = "hybrid_interface_softening_sweep_manifest_v1"
SCHEMA_SLURM_ROUND = "hybrid_interface_softening_sweep_slurm_round_v1"
SCHEMA_TASK_RESULT = "hybrid_interface_softening_sweep_task_result_v1"
SCHEMA_ANALYSIS_MANIFEST = "hybrid_interface_softening_sweep_analysis_manifest_v1"
SCHEMA_ANALYSIS_SLURM = "hybrid_interface_softening_sweep_analysis_slurm_v1"
SCHEMA_ANALYSIS_RESULT = "hybrid_interface_softening_sweep_analysis_result_v1"

DEFAULT_PDB_ID = "bilayer"
DEFAULT_LJ_ALPHAS = [0.0, 0.025, 0.05, 0.10, 0.20]
DEFAULT_SLATER_ALPHAS = [0.0, 0.25, 0.50, 1.00, 2.00]
DEFAULT_REPLICATES = 3
DEFAULT_SEED = 20260413
DEFAULT_INTEGRATION_PS_PER_STEP = 40.0
BAR_1_TO_EUP_PER_A3 = "0.000020659477"
COMP_3E4_BAR_INV_TO_A3_PER_EUP = "14.521180763676"

DEFAULT_RUNTIME_ENV = {
    "TEMPERATURE": "0.8647",
    "THERMOSTAT_TIMESCALE": "5.0",
    "THERMOSTAT_INTERVAL": "-1",
    "MIN_60_MAX_ITER": "500",
    "MIN_61_MAX_ITER": "500",
    "EQ_62_NSTEPS": "500",
    "EQ_63_NSTEPS": "500",
    "EQ_64_NSTEPS": "500",
    "EQ_65_NSTEPS": "500",
    "EQ_66_NSTEPS": "500",
    "PROD_70_NSTEPS": "10000",
    "MIN_TIME_STEP": "0.010",
    "EQ_TIME_STEP": "0.010",
    "PROD_TIME_STEP": "0.010",
    "EQ_FRAME_STEPS": "1000",
    "PROD_FRAME_STEPS": "50",
    "PROD_70_NPT_ENABLE": "1",
    "PROD_70_BAROSTAT_TYPE": "1",
    "SALT_MOLAR": "0.15",
    "ION_CUTOFF": "4.0",
    "BOX_PADDING_XY": "0.0",
    "BOX_PADDING_Z": "20.0",
    "UPSIDE_MARTINI_ENERGY_CONVERSION": "2.914952774272",
    "UPSIDE_MARTINI_LENGTH_CONVERSION": "10",
    "UPSIDE_EWALD_ENABLE": "1",
    "UPSIDE_EWALD_ALPHA": "0.2",
    "UPSIDE_EWALD_KMAX": "5",
    "UPSIDE_NPT_TAU": "4.0",
    "UPSIDE_NPT_COMPRESSIBILITY": COMP_3E4_BAR_INV_TO_A3_PER_EUP,
    "UPSIDE_NPT_COMPRESSIBILITY_XY": COMP_3E4_BAR_INV_TO_A3_PER_EUP,
    "UPSIDE_NPT_COMPRESSIBILITY_Z": "0.0",
    "UPSIDE_NPT_INTERVAL": "10",
    "UPSIDE_NPT_SEMI": "1",
    "UPSIDE_NPT_DEBUG": "1",
    "BILAYER_PDB": str(REPO_ROOT / "parameters" / "dryMARTINI" / "DOPC.pdb"),
}

DEFAULT_RUNTIME_ENV_KEYS = [
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
    "MIN_TIME_STEP",
    "EQ_TIME_STEP",
    "PROD_TIME_STEP",
    "EQ_FRAME_STEPS",
    "PROD_FRAME_STEPS",
    "PROD_70_NPT_ENABLE",
    "PROD_70_BAROSTAT_TYPE",
    "BILAYER_PDB",
    "SALT_MOLAR",
    "ION_CUTOFF",
    "BOX_PADDING_XY",
    "BOX_PADDING_Z",
    "PREP_SEED",
    "UPSIDE_MARTINI_ENERGY_CONVERSION",
    "UPSIDE_MARTINI_LENGTH_CONVERSION",
    "UPSIDE_EWALD_ENABLE",
    "UPSIDE_EWALD_ALPHA",
    "UPSIDE_EWALD_KMAX",
    "UPSIDE_NPT_TAU",
    "UPSIDE_NPT_COMPRESSIBILITY",
    "UPSIDE_NPT_COMPRESSIBILITY_XY",
    "UPSIDE_NPT_COMPRESSIBILITY_Z",
    "UPSIDE_NPT_INTERVAL",
    "UPSIDE_NPT_SEMI",
    "UPSIDE_NPT_DEBUG",
]

STAGE_SPECS = [
    {
        "label": "6.0",
        "prepare_stage": "minimization",
        "stage_kind": "minimization",
        "npt_enable": 1,
        "barostat_type": 0,
        "lipidhead_fc": 0.0,
        "nsteps_key": "MIN_60_MAX_ITER",
    },
    {
        "label": "6.1",
        "prepare_stage": "npt_prod",
        "stage_kind": "minimization",
        "npt_enable": 1,
        "barostat_type": 0,
        "lipidhead_fc": 0.0,
        "nsteps_key": "MIN_61_MAX_ITER",
    },
    {
        "label": "6.2",
        "prepare_stage": "npt_equil",
        "stage_kind": "md",
        "npt_enable": 1,
        "barostat_type": 0,
        "lipidhead_fc": 200.0,
        "nsteps_key": "EQ_62_NSTEPS",
        "dt_key": "EQ_TIME_STEP",
        "frame_key": "EQ_FRAME_STEPS",
    },
    {
        "label": "6.3",
        "prepare_stage": "npt_equil_reduced",
        "stage_kind": "md",
        "npt_enable": 1,
        "barostat_type": 0,
        "lipidhead_fc": 100.0,
        "nsteps_key": "EQ_63_NSTEPS",
        "dt_key": "EQ_TIME_STEP",
        "frame_key": "EQ_FRAME_STEPS",
    },
    {
        "label": "6.4",
        "prepare_stage": "npt_prod",
        "stage_kind": "md",
        "npt_enable": 1,
        "barostat_type": 0,
        "lipidhead_fc": 50.0,
        "nsteps_key": "EQ_64_NSTEPS",
        "dt_key": "EQ_TIME_STEP",
        "frame_key": "EQ_FRAME_STEPS",
    },
    {
        "label": "6.5",
        "prepare_stage": "npt_prod",
        "stage_kind": "md",
        "npt_enable": 1,
        "barostat_type": 0,
        "lipidhead_fc": 20.0,
        "nsteps_key": "EQ_65_NSTEPS",
        "dt_key": "EQ_TIME_STEP",
        "frame_key": "EQ_FRAME_STEPS",
    },
    {
        "label": "6.6",
        "prepare_stage": "npt_prod",
        "stage_kind": "md",
        "npt_enable": 1,
        "barostat_type": 0,
        "lipidhead_fc": 10.0,
        "nsteps_key": "EQ_66_NSTEPS",
        "dt_key": "EQ_TIME_STEP",
        "frame_key": "EQ_FRAME_STEPS",
    },
    {
        "label": "7.0",
        "prepare_stage": "npt_prod",
        "stage_kind": "md",
        "npt_enable_env_key": "PROD_70_NPT_ENABLE",
        "barostat_type_env_key": "PROD_70_BAROSTAT_TYPE",
        "lipidhead_fc": 0.0,
        "nsteps_key": "PROD_70_NSTEPS",
        "dt_key": "PROD_TIME_STEP",
        "frame_key": "PROD_FRAME_STEPS",
        "apply_softening": True,
    },
]


@dataclass
class Config:
    base_dir: Path
    pdb_id: str
    lj_alphas: List[float]
    slater_alphas: List[float]
    replicates: int
    seed: int
    integration_ps_per_step: float
    target_diffusion_um2_s: float | None
    runtime_env: Dict[str, str]


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


def _analysis_dir(base_dir: Path) -> Path:
    return base_dir / "analysis"


def _analysis_manifest_path(base_dir: Path) -> Path:
    return _analysis_dir(base_dir) / "analysis_manifest.json"


def _analysis_results_dir(base_dir: Path) -> Path:
    return _analysis_dir(base_dir) / "results" / "tasks"


def _analysis_result_path(base_dir: Path, code: str) -> Path:
    return _analysis_results_dir(base_dir) / f"{code}.json"


def _analysis_assembled_dir(base_dir: Path) -> Path:
    return _analysis_dir(base_dir) / "assembled"


def _analysis_slurm_dir(base_dir: Path) -> Path:
    return _analysis_dir(base_dir) / "slurm"


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


def _load_analysis_manifest(base_dir: Path) -> Dict[str, Any]:
    path = _analysis_manifest_path(base_dir)
    if not path.exists():
        raise RuntimeError(f"Analysis manifest not found: {path}. Run init-analysis first.")
    manifest = _load_json(path)
    if manifest.get("schema") != SCHEMA_ANALYSIS_MANIFEST:
        raise RuntimeError(f"Unexpected analysis manifest schema in {path}: {manifest.get('schema')}")
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


def _floatless_csv(text: str) -> List[str]:
    out: List[str] = []
    for chunk in text.split(","):
        item = chunk.strip()
        if item:
            out.append(item)
    return out


def _normalize_env_key_list(values: Iterable[str]) -> List[str]:
    out: List[str] = []
    seen: set[str] = set()
    for raw in values:
        key = raw.strip()
        if not key or key in seen:
            continue
        seen.add(key)
        out.append(key)
    return out


def _resolved_runtime_env() -> Dict[str, str]:
    extra = _floatless_csv(os.environ.get("HYBRID_SWEEP_EXTRA_ENV_KEYS", ""))
    keys = _normalize_env_key_list(DEFAULT_RUNTIME_ENV_KEYS + extra)
    resolved = dict(DEFAULT_RUNTIME_ENV)
    for key in keys:
        value = os.environ.get(key)
        if value is None or value == "":
            continue
        resolved[key] = value
    for key in extra:
        value = os.environ.get(key)
        if value is None or value == "":
            continue
        resolved[key] = value
    return resolved


def _build_config(args: argparse.Namespace, base_dir: Path) -> Config:
    lj_alphas = sorted(dict.fromkeys(float(x) for x in args.lj_alphas))
    slater_alphas = sorted(dict.fromkeys(float(x) for x in args.slater_alphas))
    for value in lj_alphas:
        if not math.isfinite(value) or value < 0.0:
            raise ValueError(f"lj alpha must be finite and >= 0, got {value}")
    for value in slater_alphas:
        if not math.isfinite(value) or value < 0.0:
            raise ValueError(f"slater alpha must be finite and >= 0, got {value}")
    replicates = int(args.replicates)
    if replicates < 1:
        raise ValueError(f"replicates must be >= 1, got {replicates}")
    integration_ps_per_step = float(args.integration_ps_per_step)
    if not math.isfinite(integration_ps_per_step) or integration_ps_per_step <= 0.0:
        raise ValueError(
            f"integration_ps_per_step must be finite and > 0, got {integration_ps_per_step}"
        )
    target = args.target_diffusion_um2_s
    if target is not None:
        target = float(target)
        if not math.isfinite(target) or target <= 0.0:
            raise ValueError(f"target diffusion must be finite and > 0, got {target}")
    return Config(
        base_dir=base_dir,
        pdb_id=str(args.pdb_id),
        lj_alphas=lj_alphas,
        slater_alphas=slater_alphas,
        replicates=replicates,
        seed=int(args.seed),
        integration_ps_per_step=integration_ps_per_step,
        target_diffusion_um2_s=target,
        runtime_env=_resolved_runtime_env(),
    )


def _build_manifest(config: Config) -> Dict[str, Any]:
    tasks: List[Dict[str, Any]] = []
    task_id = 0
    for lj_alpha in config.lj_alphas:
        for slater_alpha in config.slater_alphas:
            for replicate in range(config.replicates):
                code = (
                    f"lj{_format_float_tag(lj_alpha)}_"
                    f"coul{_format_float_tag(slater_alpha)}_"
                    f"r{replicate + 1:02d}"
                )
                tasks.append(
                    {
                        "task_id": task_id,
                        "code": code,
                        "lj_alpha": float(lj_alpha),
                        "slater_alpha": float(slater_alpha),
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
            "lj_alphas": config.lj_alphas,
            "slater_alphas": config.slater_alphas,
            "replicates": config.replicates,
            "seed": config.seed,
            "prep_script": str(PREP_SCRIPT),
            "upside_executable": str(UPSIDE_EXECUTABLE),
            "runtime_env": config.runtime_env,
            "integration_ps_per_step": config.integration_ps_per_step,
            "target_diffusion_um2_s": config.target_diffusion_um2_s,
        },
        "tasks": tasks,
    }


def _require_runtime_files(runtime_env: Dict[str, str]) -> None:
    if not PREP_SCRIPT.exists():
        raise RuntimeError(f"Preparation script not found: {PREP_SCRIPT}")
    if not UPSIDE_EXECUTABLE.exists():
        raise RuntimeError(f"Upside executable not found: {UPSIDE_EXECUTABLE}")
    bilayer_pdb = Path(runtime_env["BILAYER_PDB"]).expanduser().resolve()
    if not bilayer_pdb.exists():
        raise RuntimeError(f"Bilayer PDB not found: {bilayer_pdb}")


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


def _stage_file_path(checkpoint_dir: Path, pdb_id: str, label: str, prepared: bool) -> Path:
    suffix = "prepared.up" if prepared else "up"
    return checkpoint_dir / f"{pdb_id}.stage_{label}.{suffix}"


def _prepare_runtime_structure(
    *,
    pdb_id: str,
    runtime_pdb: Path,
    runtime_itp: Path,
    run_dir: Path,
    summary_json: Path,
    env: Dict[str, str],
    log_path: Path,
) -> None:
    cmd = [
        "python3",
        str(PREP_SCRIPT),
        "--mode",
        "bilayer",
        "--pdb-id",
        pdb_id,
        "--runtime-pdb-output",
        str(runtime_pdb),
        "--runtime-itp-output",
        str(runtime_itp),
        "--prepare-structure",
        "1",
        "--run-dir",
        str(run_dir),
        "--bilayer-pdb",
        env["BILAYER_PDB"],
        "--salt-molar",
        env["SALT_MOLAR"],
        "--ion-cutoff",
        env["ION_CUTOFF"],
        "--box-padding-xy",
        env["BOX_PADDING_XY"],
        "--box-padding-z",
        env["BOX_PADDING_Z"],
        "--seed",
        env["PREP_SEED"],
        "--summary-json",
        str(summary_json),
    ]
    _call_logged(cmd, log_path, env)


def _stage_env(env: Dict[str, str], stage_label: str, npt_enable: int, lipidhead_fc: float) -> Dict[str, str]:
    out = dict(env)
    out["UPSIDE_NPT_ENABLE"] = str(int(npt_enable))
    out["UPSIDE_BILAYER_LIPIDHEAD_FC"] = f"{float(lipidhead_fc):.8g}"
    if not npt_enable:
        return out
    if stage_label in {"6.0", "6.1"}:
        target = "0.0"
    else:
        target = BAR_1_TO_EUP_PER_A3
    out["UPSIDE_NPT_TARGET_PXY"] = target
    out["UPSIDE_NPT_TARGET_PZ"] = target
    return out


def _prepare_stage_input(
    *,
    pdb_id: str,
    runtime_pdb: Path,
    runtime_itp: Path,
    run_dir: Path,
    prepare_stage: str,
    summary_json: Path,
    env: Dict[str, str],
    log_path: Path,
) -> None:
    stage_env = dict(env)
    stage_env["UPSIDE_SIMULATION_STAGE"] = prepare_stage
    cmd = [
        "python3",
        str(PREP_SCRIPT),
        "--mode",
        "bilayer",
        "--pdb-id",
        pdb_id,
        "--runtime-pdb-output",
        str(runtime_pdb),
        "--runtime-itp-output",
        str(runtime_itp),
        "--prepare-structure",
        "0",
        "--stage",
        prepare_stage,
        "--run-dir",
        str(run_dir),
        "--summary-json",
        str(summary_json),
    ]
    _call_logged(cmd, log_path, stage_env)


def _copy_last_position(input_file: Path, output_file: Path) -> None:
    with h5py.File(input_file, "r") as src:
        if "/output/pos" in src and src["/output/pos"].shape[0] > 0:
            last_pos = np.asarray(src["/output/pos"][-1, 0, :, :], dtype=np.float64)
            last_pos = last_pos[:, :, np.newaxis]
        else:
            last_pos = np.asarray(src["/input/pos"][:, :, 0], dtype=np.float64)
            last_pos = last_pos[:, :, np.newaxis]

        last_box = None
        if "/output/box" in src:
            box_data = np.asarray(src["/output/box"][:], dtype=np.float64)
            if box_data.size > 0:
                box_frame = box_data[-1]
                if box_frame.ndim == 2 and box_frame.shape[0] > 0:
                    last_box = np.asarray(box_frame[0], dtype=np.float64)
                elif box_frame.ndim == 1 and box_frame.shape[0] == 3:
                    last_box = np.asarray(box_frame, dtype=np.float64)
        if last_box is None and "/input/potential/martini_potential" in src:
            pot = src["/input/potential/martini_potential"]
            if all(key in pot.attrs for key in ("x_len", "y_len", "z_len")):
                last_box = np.asarray(
                    [pot.attrs["x_len"], pot.attrs["y_len"], pot.attrs["z_len"]],
                    dtype=np.float64,
                )

    with h5py.File(output_file, "r+") as dst:
        target_pos = np.asarray(dst["/input/pos"][:], dtype=np.float64)
        if target_pos.shape[0] != last_pos.shape[0]:
            merged = target_pos.copy()
            merged[: min(target_pos.shape[0], last_pos.shape[0]), :, :] = last_pos[
                : min(target_pos.shape[0], last_pos.shape[0]),
                :,
                :,
            ]
            last_pos = merged
        if "/input/pos" in dst:
            del dst["/input/pos"]
        dst.create_dataset("/input/pos", data=last_pos)
        if last_box is not None and "/input/potential/martini_potential" in dst:
            pot = dst["/input/potential/martini_potential"]
            pot.attrs["x_len"] = float(last_box[0])
            pot.attrs["y_len"] = float(last_box[1])
            pot.attrs["z_len"] = float(last_box[2])


def _set_barostat_type(up_file: Path, barostat_type: int) -> None:
    with h5py.File(up_file, "r+") as h5:
        if "/input/barostat" not in h5:
            return
        h5["/input/barostat"].attrs["type"] = int(barostat_type)


def _set_production_softening(up_file: Path, *, lj_alpha: float, slater_alpha: float) -> Dict[str, Any]:
    lj_enable = int(lj_alpha > 0.0)
    coul_enable = int(slater_alpha > 0.0)
    with h5py.File(up_file, "r+") as h5:
        grp = h5["/input/potential/martini_potential"]
        grp.attrs["lj_soften"] = np.int8(lj_enable)
        grp.attrs["lj_soften_alpha"] = np.float32(float(max(0.0, lj_alpha)))
        grp.attrs["coulomb_soften"] = np.int8(coul_enable)
        grp.attrs["slater_alpha"] = np.float32(float(max(0.0, slater_alpha)))
    return {
        "lj_soften": lj_enable,
        "lj_soften_alpha": float(max(0.0, lj_alpha)),
        "coulomb_soften": coul_enable,
        "slater_alpha": float(max(0.0, slater_alpha)),
    }


def _run_minimization_stage(
    *,
    stage_label: str,
    up_file: Path,
    max_iter: int,
    env: Dict[str, str],
    log_path: Path,
) -> None:
    cmd = [
        str(UPSIDE_EXECUTABLE),
        str(up_file),
        "--duration",
        "0",
        "--frame-interval",
        "1",
        "--temperature",
        env["TEMPERATURE"],
        "--time-step",
        env["MIN_TIME_STEP"],
        "--thermostat-timescale",
        env["THERMOSTAT_TIMESCALE"],
        "--thermostat-interval",
        env["THERMOSTAT_INTERVAL"],
        "--seed",
        env["SEED"],
        "--integrator",
        "v",
        "--disable-recentering",
        "--minimize",
        "--min-max-iter",
        str(int(max_iter)),
        "--min-energy-tol",
        "1e-6",
        "--min-force-tol",
        "1e-3",
        "--min-step",
        "0.01",
    ]
    _call_logged(cmd, log_path, env)


def _run_md_stage(
    *,
    up_file: Path,
    nsteps: int,
    dt: float,
    frame_steps: int,
    env: Dict[str, str],
    log_path: Path,
) -> Dict[str, Any]:
    effective_frame_steps = int(frame_steps)
    if effective_frame_steps >= int(nsteps):
        effective_frame_steps = max(1, int(nsteps) // 10)
    frame_interval = effective_frame_steps * float(dt)
    cmd = [
        str(UPSIDE_EXECUTABLE),
        str(up_file),
        "--duration-steps",
        str(int(nsteps)),
        "--frame-interval",
        f"{frame_interval:.10g}",
        "--temperature",
        env["TEMPERATURE"],
        "--time-step",
        f"{float(dt):.10g}",
        "--thermostat-timescale",
        env["THERMOSTAT_TIMESCALE"],
        "--thermostat-interval",
        env["THERMOSTAT_INTERVAL"],
        "--seed",
        env["SEED"],
        "--integrator",
        "v",
        "--disable-recentering",
    ]
    _call_logged(cmd, log_path, env)
    return {
        "nsteps": int(nsteps),
        "time_step": float(dt),
        "frame_steps": int(effective_frame_steps),
        "frame_interval_time": float(frame_interval),
    }


def cmd_init_run(args: argparse.Namespace) -> int:
    base_dir = Path(args.base_dir).expanduser().resolve()
    base_dir.mkdir(parents=True, exist_ok=True)
    config = _build_config(args, base_dir)
    _require_runtime_files(config.runtime_env)
    manifest = _build_manifest(config)
    _write_json(_manifest_path(base_dir), manifest)
    print(f"Initialized hybrid interface softening sweep: {base_dir}")
    print(f"Task count: {len(manifest['tasks'])}")
    print(f"LJ alphas: {manifest['settings']['lj_alphas']}")
    print(f"Slater alphas: {manifest['settings']['slater_alphas']}")
    return 0


def _run_task(base_dir: Path, manifest: Dict[str, Any], task: Dict[str, Any], overwrite: bool) -> Dict[str, Any]:
    result_path = _result_path(base_dir, task["code"])
    if result_path.exists() and not overwrite:
        existing = _load_json(result_path)
        if existing.get("success"):
            print(f"Skipping completed task {task['code']}")
            return existing

    settings = manifest["settings"]
    runtime_env = dict(settings["runtime_env"])
    task_dir = _task_dir(base_dir, task["code"])
    run_dir = task_dir / "run"
    checkpoint_dir = run_dir / "checkpoints"
    log_dir = run_dir / "logs"
    input_dir = task_dir / "inputs"
    prep_dir = task_dir / "prep"
    prep_dir.mkdir(parents=True, exist_ok=True)
    runtime_pdb = input_dir / f"{settings['pdb_id']}.MARTINI.pdb"
    runtime_itp = input_dir / f"{settings['pdb_id']}_unused.itp"

    env = _base_runtime_env()
    env.update({str(k): str(v) for k, v in runtime_env.items()})
    env["SEED"] = str(int(task["seed"]))
    env.setdefault("PREP_SEED", env["SEED"])

    stage_logs: Dict[str, str] = {}
    md_runtime: Dict[str, Dict[str, Any]] = {}
    stage_70_softening: Dict[str, Any] | None = None

    try:
        _prepare_runtime_structure(
            pdb_id=str(settings["pdb_id"]),
            runtime_pdb=runtime_pdb,
            runtime_itp=runtime_itp,
            run_dir=run_dir,
            summary_json=prep_dir / "runtime_structure.summary.json",
            env=env,
            log_path=log_dir / "prepare_runtime_structure.log",
        )
        stage_logs["prepare_runtime_structure"] = str(log_dir / "prepare_runtime_structure.log")

        previous_stage_file: Path | None = None
        for spec in STAGE_SPECS:
            label = str(spec["label"])
            prepare_stage = str(spec["prepare_stage"])
            prepared_file = _stage_file_path(checkpoint_dir, str(settings["pdb_id"]), label, prepared=True)
            stage_file = _stage_file_path(checkpoint_dir, str(settings["pdb_id"]), label, prepared=False)
            npt_enable = int(runtime_env.get(spec.get("npt_enable_env_key", ""), spec.get("npt_enable", 0)))
            barostat_type = int(
                runtime_env.get(spec.get("barostat_type_env_key", ""), spec.get("barostat_type", 0))
            )
            stage_env = _stage_env(env, label, npt_enable, float(spec["lipidhead_fc"]))
            prepare_log = log_dir / f"stage_{label}.prepare.log"
            _prepare_stage_input(
                pdb_id=str(settings["pdb_id"]),
                runtime_pdb=runtime_pdb,
                runtime_itp=runtime_itp,
                run_dir=run_dir,
                prepare_stage=prepare_stage,
                summary_json=prep_dir / f"stage_{label}.summary.json",
                env=stage_env,
                log_path=prepare_log,
            )
            stage_logs[f"stage_{label}_prepare"] = str(prepare_log)

            prepared_tmp = run_dir / "test.input.up"
            if not prepared_tmp.exists():
                raise RuntimeError(f"Preparation failed for stage {label}: {prepared_tmp} not found")
            prepared_file.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(prepared_tmp), str(prepared_file))
            shutil.copy2(prepared_file, stage_file)

            if previous_stage_file is not None:
                _copy_last_position(previous_stage_file, stage_file)

            if npt_enable:
                _set_barostat_type(stage_file, barostat_type)

            if spec.get("apply_softening"):
                stage_70_softening = _set_production_softening(
                    stage_file,
                    lj_alpha=float(task["lj_alpha"]),
                    slater_alpha=float(task["slater_alpha"]),
                )

            run_log = log_dir / f"stage_{label}.run.log"
            if spec["stage_kind"] == "minimization":
                _run_minimization_stage(
                    stage_label=label,
                    up_file=stage_file,
                    max_iter=int(runtime_env[spec["nsteps_key"]]),
                    env=stage_env,
                    log_path=run_log,
                )
            else:
                md_runtime[label] = _run_md_stage(
                    up_file=stage_file,
                    nsteps=int(runtime_env[spec["nsteps_key"]]),
                    dt=float(runtime_env[spec["dt_key"]]),
                    frame_steps=int(runtime_env[spec["frame_key"]]),
                    env=stage_env,
                    log_path=run_log,
                )
            stage_logs[f"stage_{label}_run"] = str(run_log)
            previous_stage_file = stage_file

        stage_70_file = _stage_file_path(checkpoint_dir, str(settings["pdb_id"]), "7.0", prepared=False)
        if not stage_70_file.exists():
            raise RuntimeError(f"Missing expected stage-7 checkpoint: {stage_70_file}")

        result = {
            "schema": SCHEMA_TASK_RESULT,
            "created_at_utc": _now_utc(),
            "base_dir": str(base_dir),
            "task": task,
            "success": True,
            "run_dir": str(run_dir),
            "checkpoint_dir": str(checkpoint_dir),
            "log_dir": str(log_dir),
            "runtime_pdb": str(runtime_pdb),
            "stage_logs": stage_logs,
            "stage_70_file": str(stage_70_file),
            "production_softening": stage_70_softening,
            "production_time_step": float(runtime_env["PROD_TIME_STEP"]),
            "integration_ps_per_step": float(settings["integration_ps_per_step"]),
            "target_diffusion_um2_s": settings["target_diffusion_um2_s"],
            "md_runtime": md_runtime,
            "runtime_env": runtime_env,
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
            "log_dir": str(log_dir),
            "error": str(exc),
            "stage_logs": stage_logs,
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
    failures: List[str] = []
    for task in selected:
        try:
            _run_task(base_dir, manifest, task, overwrite=bool(args.overwrite))
        except Exception as exc:
            failures.append(f"{task['code']}: {exc}")
    if not args.no_assemble:
        assemble_results(base_dir)
    if failures:
        raise RuntimeError("Sweep failures:\n" + "\n".join(failures))
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
    grouped: Dict[tuple[float, float], List[Dict[str, Any]]] = {}
    for result in results:
        task = result["task"]
        row = {
            "task_id": int(task["task_id"]),
            "code": str(task["code"]),
            "lj_alpha": float(task["lj_alpha"]),
            "slater_alpha": float(task["slater_alpha"]),
            "replicate": int(task["replicate"]),
            "seed": int(task["seed"]),
            "production_time_step": float(result["production_time_step"]),
            "integration_ps_per_step": float(result["integration_ps_per_step"]),
            "run_dir": str(result["run_dir"]),
            "log_dir": str(result["log_dir"]),
            "stage_70_file": str(result["stage_70_file"]),
        }
        task_rows.append(row)
        grouped.setdefault((row["lj_alpha"], row["slater_alpha"]), []).append(row)

    task_csv = assembled_dir / "task_results.csv"
    with task_csv.open("w", encoding="utf-8", newline="") as fh:
        fieldnames = [
            "task_id",
            "code",
            "lj_alpha",
            "slater_alpha",
            "replicate",
            "seed",
            "production_time_step",
            "integration_ps_per_step",
            "run_dir",
            "log_dir",
            "stage_70_file",
        ]
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in sorted(task_rows, key=lambda item: item["task_id"]):
            writer.writerow(row)

    summary_rows: List[Dict[str, Any]] = []
    expected_replicates = int(manifest["settings"]["replicates"])
    for (lj_alpha, slater_alpha), rows in sorted(grouped.items()):
        summary_rows.append(
            {
                "lj_alpha": float(lj_alpha),
                "slater_alpha": float(slater_alpha),
                "n_replicates_expected": expected_replicates,
                "n_replicates_completed": len(rows),
                "all_stage70_present": int(all(Path(row["stage_70_file"]).exists() for row in rows)),
            }
        )

    summary_csv = assembled_dir / "condition_summary.csv"
    with summary_csv.open("w", encoding="utf-8", newline="") as fh:
        fieldnames = [
            "lj_alpha",
            "slater_alpha",
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
        "n_conditions_completed": len(summary_rows),
        "task_results_csv": str(task_csv),
        "condition_summary_csv": str(summary_csv),
    }
    _write_json(assembled_dir / "summary.json", payload)
    print(f"Assembled sweep results written under {assembled_dir}")
    return 0


def cmd_assemble_results(args: argparse.Namespace) -> int:
    base_dir = Path(args.base_dir).expanduser().resolve()
    return assemble_results(base_dir)


def _normalize_output_positions(arr: np.ndarray) -> np.ndarray:
    if arr.ndim == 4:
        if arr.shape[-1] == 3:
            return np.asarray(arr[:, 0, :, :], dtype=np.float64)
        if arr.shape[1] == 3:
            return np.asarray(np.transpose(arr[:, :, :, 0], (0, 2, 1)), dtype=np.float64)
    if arr.ndim == 3:
        if arr.shape[-1] == 3:
            return np.asarray(arr, dtype=np.float64)
        if arr.shape[1] == 3:
            return np.asarray(np.transpose(arr, (0, 2, 1)), dtype=np.float64)
    raise ValueError(f"Unexpected output/pos shape: {arr.shape}")


def _normalize_time(arr: np.ndarray) -> np.ndarray:
    out = np.asarray(arr)
    while out.ndim > 1:
        out = out[:, 0]
    return np.asarray(out, dtype=np.float64).reshape(-1)


def _normalize_box(arr: np.ndarray) -> np.ndarray:
    out = np.asarray(arr)
    if out.ndim == 3:
        if out.shape[-1] == 3:
            return np.asarray(out[:, 0, :], dtype=np.float64)
        if out.shape[1] == 3:
            return np.asarray(np.transpose(out[:, :, 0], (0, 1)), dtype=np.float64)
    if out.ndim == 2 and out.shape[-1] == 3:
        return np.asarray(out, dtype=np.float64)
    raise ValueError(f"Unexpected output/box shape: {arr.shape}")


def _decode_str_array(dataset: h5py.Dataset) -> np.ndarray:
    arr = dataset[:]
    out = []
    for item in arr:
        if isinstance(item, (bytes, np.bytes_)):
            out.append(item.decode("utf-8", errors="ignore").strip())
        else:
            out.append(str(item).strip())
    return np.asarray(out, dtype=object)


def _load_box_frames(h5: h5py.File, n_frame: int) -> np.ndarray:
    if "output/box" in h5:
        return _normalize_box(h5["output/box"][:])
    grp = h5["input/potential/martini_potential"]
    box = np.array([grp.attrs["x_len"], grp.attrs["y_len"], grp.attrs["z_len"]], dtype=np.float64)
    return np.repeat(box[None, :], n_frame, axis=0)


def _minimum_image_xy(delta_xy: np.ndarray, box_xy: np.ndarray) -> np.ndarray:
    return delta_xy - box_xy * np.round(delta_xy / box_xy)


def _unwrap_xy(positions_xy: np.ndarray, box_xy: np.ndarray) -> np.ndarray:
    unwrapped = np.zeros_like(positions_xy, dtype=np.float64)
    unwrapped[0] = positions_xy[0]
    for frame in range(1, positions_xy.shape[0]):
        delta = positions_xy[frame] - positions_xy[frame - 1]
        delta = _minimum_image_xy(delta, box_xy[frame - 1][None, :])
        unwrapped[frame] = unwrapped[frame - 1] + delta
    return unwrapped


def _compute_msd(xy: np.ndarray, times: np.ndarray, max_lag: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    msd = np.zeros((max_lag,), dtype=np.float64)
    lag_times = np.zeros((max_lag,), dtype=np.float64)
    counts = np.zeros((max_lag,), dtype=np.int64)
    for lag in range(1, max_lag + 1):
        disp = xy[lag:] - xy[:-lag]
        sq = np.sum(disp * disp, axis=2)
        msd[lag - 1] = float(np.mean(sq))
        lag_times[lag - 1] = float(np.mean(times[lag:] - times[:-lag]))
        counts[lag - 1] = int(sq.size)
    return lag_times, msd, counts


def _fit_diffusion(lag_times: np.ndarray, msd: np.ndarray) -> Dict[str, Any]:
    max_lag = lag_times.shape[0]
    if max_lag < 5:
        raise ValueError("Not enough lag points to fit diffusion")
    fit_start = max(1, int(math.ceil(0.10 * max_lag)))
    fit_end = max(fit_start + 2, int(math.floor(0.40 * max_lag)))
    fit_slice = slice(fit_start - 1, fit_end)
    x = lag_times[fit_slice]
    y = msd[fit_slice]
    if x.shape[0] < 2:
        raise ValueError("Not enough points in diffusion fit window")
    slope, intercept = np.polyfit(x, y, 1)
    y_pred = slope * x + intercept
    ss_res = float(np.sum((y - y_pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 if ss_tot <= 0.0 else 1.0 - ss_res / ss_tot
    return {
        "fit_start_lag_index_1based": fit_start,
        "fit_end_lag_index_1based": fit_end,
        "slope_angstrom2_per_time": float(slope),
        "intercept_angstrom2": float(intercept),
        "r2": float(r2),
        "n_fit_points": int(x.shape[0]),
        "diffusion_angstrom2_per_time": float(slope / 4.0),
        "diffusion_nm2_per_time": float(slope / 400.0),
    }


def _select_lipid_atoms_and_po4(
    atom_names: np.ndarray,
    residue_names: np.ndarray | None,
    molecule_ids: np.ndarray | None,
    particle_class: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray]:
    if particle_class is not None:
        lipid_mask = particle_class == "LIPID"
        po4_mask = lipid_mask & (atom_names == "PO4")
        if np.any(po4_mask):
            return lipid_mask, po4_mask

    if residue_names is not None:
        lipid_mask = np.isin(residue_names, np.asarray(["DOP", "DOPC"], dtype=object))
        po4_mask = lipid_mask & (atom_names == "PO4")
        if np.any(po4_mask):
            return lipid_mask, po4_mask

    po4_mask = atom_names == "PO4"
    if molecule_ids is None or not np.any(po4_mask):
        raise ValueError("Could not identify bilayer PO4 beads for analysis")
    lipid_molecule_ids = np.unique(molecule_ids[po4_mask])
    lipid_mask = np.isin(molecule_ids, lipid_molecule_ids)
    po4_mask = lipid_mask & (atom_names == "PO4")
    if not np.any(po4_mask):
        raise ValueError("Could not infer lipid molecules from PO4 beads")
    return lipid_mask, po4_mask


def _convert_diffusion_units(
    diffusion_nm2_per_time: float,
    production_time_step: float,
    integration_ps_per_step: float,
) -> Dict[str, float]:
    if not (production_time_step > 0.0 and integration_ps_per_step > 0.0):
        raise ValueError("Production timestep and integration ps/step must be > 0")
    ns_per_time_unit = integration_ps_per_step / production_time_step / 1000.0
    diffusion_nm2_per_ns = diffusion_nm2_per_time / ns_per_time_unit
    diffusion_um2_per_s = diffusion_nm2_per_ns * 1000.0
    viscosity_proxy_s_per_um2 = float("inf")
    if diffusion_um2_per_s > 0.0:
        viscosity_proxy_s_per_um2 = 1.0 / diffusion_um2_per_s
    return {
        "ns_per_time_unit": float(ns_per_time_unit),
        "diffusion_nm2_per_ns": float(diffusion_nm2_per_ns),
        "diffusion_um2_per_s": float(diffusion_um2_per_s),
        "viscosity_proxy_s_per_um2": float(viscosity_proxy_s_per_um2),
    }


def _analyze_stage7(
    stage_file: Path,
    *,
    production_time_step: float,
    integration_ps_per_step: float,
    target_diffusion_um2_s: float | None,
) -> Dict[str, Any]:
    with h5py.File(stage_file, "r") as h5:
        if "output/pos" not in h5:
            raise ValueError(f"Missing output/pos in {stage_file}")
        pos = _normalize_output_positions(h5["output/pos"][:])
        n_frame, n_atom, _ = pos.shape
        if n_frame < 10:
            raise ValueError(f"Need at least 10 output frames for MSD analysis, found {n_frame}")

        if "output/time" in h5:
            times = _normalize_time(h5["output/time"][:])
        else:
            times = np.arange(n_frame, dtype=np.float64)
        if times.shape[0] != n_frame:
            raise ValueError("output/time length does not match output/pos frames")

        box = _load_box_frames(h5, n_frame)
        inp = h5["input"]
        atom_names = _decode_str_array(inp["atom_names"])
        residue_names = _decode_str_array(inp["residue_names"]) if "residue_names" in inp else None
        molecule_ids = np.asarray(inp["molecule_ids"][:], dtype=np.int32) if "molecule_ids" in inp else None
        particle_class = _decode_str_array(inp["particle_class"]) if "particle_class" in inp else None

        lipid_mask, po4_mask = _select_lipid_atoms_and_po4(atom_names, residue_names, molecule_ids, particle_class)
        lipid_indices = np.where(lipid_mask)[0]
        po4_indices = np.where(po4_mask)[0]
        if lipid_indices.size == 0 or po4_indices.size == 0:
            raise ValueError("No lipid PO4 atoms selected for bilayer analysis")

        burn_start = max(1, int(math.floor(0.20 * n_frame)))
        pos = pos[burn_start:]
        times = times[burn_start:]
        box = box[burn_start:]
        n_frame = pos.shape[0]
        if n_frame < 10:
            raise ValueError("Too few frames remain after burn-in for MSD analysis")

        box_xy = box[:, :2]
        lipid_xy = pos[:, lipid_indices, :2]
        po4_xy = pos[:, po4_indices, :2]

        lipid_xy_unwrapped = _unwrap_xy(lipid_xy, box_xy)
        po4_xy_unwrapped = _unwrap_xy(po4_xy, box_xy)
        bilayer_com_xy = np.mean(lipid_xy_unwrapped, axis=1)
        po4_xy_corrected = po4_xy_unwrapped - bilayer_com_xy[:, None, :]

        max_lag = min(max(5, int(math.floor(0.25 * n_frame))), n_frame - 1)
        if max_lag < 5:
            raise ValueError("Not enough post-burn-in frames for lag-window analysis")

        lag_times, msd_po4, counts_po4 = _compute_msd(po4_xy_corrected, times, max_lag)
        fit_po4 = _fit_diffusion(lag_times, msd_po4)
        converted = _convert_diffusion_units(
            fit_po4["diffusion_nm2_per_time"],
            production_time_step=production_time_step,
            integration_ps_per_step=integration_ps_per_step,
        )

        analysis = {
            "n_frames_total": int(burn_start + n_frame),
            "n_frames_used": int(n_frame),
            "burn_in_frames": int(burn_start),
            "n_atom": int(n_atom),
            "n_lipid_atoms": int(lipid_indices.shape[0]),
            "n_po4": int(po4_indices.shape[0]),
            "msd_lag_time": lag_times.tolist(),
            "po4_msd_angstrom2": msd_po4.tolist(),
            "po4_msd_counts": counts_po4.tolist(),
            "fit_po4": fit_po4,
            "po4_lateral_diffusion_angstrom2_per_time": float(fit_po4["diffusion_angstrom2_per_time"]),
            "po4_lateral_diffusion_nm2_per_time": float(fit_po4["diffusion_nm2_per_time"]),
            "po4_lateral_diffusion_nm2_per_ns": float(converted["diffusion_nm2_per_ns"]),
            "po4_lateral_diffusion_um2_per_s": float(converted["diffusion_um2_per_s"]),
            "viscosity_proxy_s_per_um2": float(converted["viscosity_proxy_s_per_um2"]),
            "integration_ps_per_step": float(integration_ps_per_step),
            "production_time_step": float(production_time_step),
            "ns_per_time_unit": float(converted["ns_per_time_unit"]),
        }
        if target_diffusion_um2_s is not None:
            analysis["target_diffusion_um2_per_s"] = float(target_diffusion_um2_s)
            analysis["target_abs_error_um2_per_s"] = float(
                abs(analysis["po4_lateral_diffusion_um2_per_s"] - float(target_diffusion_um2_s))
            )
        return analysis


def _discover_analysis_tasks(base_dir: Path, manifest: Dict[str, Any]) -> List[Dict[str, Any]]:
    task_by_code = {str(task["code"]): dict(task) for task in manifest["tasks"]}
    discovered: List[Dict[str, Any]] = []
    for result in _successful_results_from_dir(_result_dir(base_dir)):
        task = dict(result["task"])
        code = str(task["code"])
        task.update(task_by_code.get(code, {}))
        discovered.append(
            {
                "analysis_task_id": len(discovered),
                "code": code,
                "task": task,
                "stage7_file": str(result["stage_70_file"]),
                "production_time_step": float(result["production_time_step"]),
                "integration_ps_per_step": float(result["integration_ps_per_step"]),
                "target_diffusion_um2_s": result.get("target_diffusion_um2_s"),
            }
        )
    if not discovered:
        raise RuntimeError(f"No successful sweep task results found under {_result_dir(base_dir)}")
    discovered.sort(key=lambda item: (int(item["task"].get("task_id", 10**9)), item["code"]))
    for idx, item in enumerate(discovered):
        item["analysis_task_id"] = idx
    return discovered


def cmd_init_analysis(args: argparse.Namespace) -> int:
    base_dir = Path(args.base_dir).expanduser().resolve()
    manifest = _load_manifest(base_dir)
    tasks = _discover_analysis_tasks(base_dir, manifest)
    payload = {
        "schema": SCHEMA_ANALYSIS_MANIFEST,
        "created_at_utc": _now_utc(),
        "base_dir": str(base_dir),
        "manifest_path": str(_manifest_path(base_dir)),
        "pdb_id": manifest["settings"]["pdb_id"],
        "n_stage7_files": len(tasks),
        "tasks": tasks,
    }
    _write_json(_analysis_manifest_path(base_dir), payload)
    print(f"Analysis manifest: {_analysis_manifest_path(base_dir)}")
    print(f"Discovered stage_7.0 files: {len(tasks)}")
    return 0


def _run_analysis_task(base_dir: Path, task_entry: Dict[str, Any], overwrite: bool) -> Dict[str, Any]:
    task = dict(task_entry["task"])
    code = str(task["code"])
    result_path = _analysis_result_path(base_dir, code)
    if result_path.exists() and not overwrite:
        existing = _load_json(result_path)
        if existing.get("success"):
            print(f"Skipping completed analysis task {code}")
            return existing

    stage7_file = Path(task_entry["stage7_file"]).expanduser().resolve()
    if not stage7_file.exists():
        raise RuntimeError(f"Stage-7 checkpoint not found: {stage7_file}")

    try:
        analysis = _analyze_stage7(
            stage7_file,
            production_time_step=float(task_entry["production_time_step"]),
            integration_ps_per_step=float(task_entry["integration_ps_per_step"]),
            target_diffusion_um2_s=(
                None
                if task_entry.get("target_diffusion_um2_s") is None
                else float(task_entry["target_diffusion_um2_s"])
            ),
        )
        result = {
            "schema": SCHEMA_ANALYSIS_RESULT,
            "created_at_utc": _now_utc(),
            "base_dir": str(base_dir),
            "analysis_task_id": int(task_entry["analysis_task_id"]),
            "success": True,
            "task": task,
            "stage7_file": str(stage7_file),
            "analysis": analysis,
        }
        _write_json(result_path, result)
        return result
    except Exception as exc:
        result = {
            "schema": SCHEMA_ANALYSIS_RESULT,
            "created_at_utc": _now_utc(),
            "base_dir": str(base_dir),
            "analysis_task_id": int(task_entry["analysis_task_id"]),
            "success": False,
            "task": task,
            "stage7_file": str(stage7_file),
            "error": str(exc),
        }
        _write_json(result_path, result)
        raise


def cmd_run_analysis_local(args: argparse.Namespace) -> int:
    base_dir = Path(args.base_dir).expanduser().resolve()
    analysis_manifest = _load_analysis_manifest(base_dir)
    tasks = analysis_manifest["tasks"]
    start_task = max(0, int(args.start_task))
    selected = tasks[start_task:]
    if args.max_tasks is not None:
        selected = selected[: max(0, int(args.max_tasks))]

    failures: List[str] = []
    for task_entry in selected:
        try:
            _run_analysis_task(base_dir, task_entry, overwrite=bool(args.overwrite))
        except Exception as exc:
            failures.append(f"{task_entry['code']}: {exc}")

    if not args.no_assemble:
        assemble_analysis(base_dir)

    if failures:
        raise RuntimeError("Analysis failures:\n" + "\n".join(failures))
    return 0


def cmd_run_analysis_task(args: argparse.Namespace) -> int:
    round_manifest = _load_json(Path(args.round_manifest).expanduser().resolve())
    if round_manifest.get("schema") != SCHEMA_ANALYSIS_SLURM:
        raise RuntimeError(f"Unexpected analysis round manifest schema: {round_manifest.get('schema')}")
    base_dir = Path(round_manifest["base_dir"]).expanduser().resolve()
    analysis_manifest = _load_analysis_manifest(base_dir)
    task_id = int(args.task_id)
    tasks = analysis_manifest["tasks"]
    if task_id < 0 or task_id >= len(tasks):
        raise RuntimeError(f"Analysis task id out of range: {task_id} for {len(tasks)} tasks")
    _run_analysis_task(base_dir, tasks[task_id], overwrite=bool(args.overwrite))
    return 0


def assemble_analysis(base_dir: Path) -> int:
    manifest = _load_manifest(base_dir)
    analysis_manifest = _load_analysis_manifest(base_dir)
    results_dir = _analysis_results_dir(base_dir)
    assembled_dir = _analysis_assembled_dir(base_dir)
    assembled_dir.mkdir(parents=True, exist_ok=True)

    successful_results = _successful_results_from_dir(results_dir)
    task_rows: List[Dict[str, Any]] = []
    grouped: Dict[tuple[float, float], List[Dict[str, Any]]] = {}
    target_diffusion_um2_s = manifest["settings"].get("target_diffusion_um2_s")
    for result in successful_results:
        task = result["task"]
        analysis = result["analysis"]
        row = {
            "analysis_task_id": int(result["analysis_task_id"]),
            "task_id": int(task.get("task_id", -1)),
            "code": str(task["code"]),
            "lj_alpha": float(task["lj_alpha"]),
            "slater_alpha": float(task["slater_alpha"]),
            "replicate": int(task["replicate"]),
            "po4_lateral_diffusion_angstrom2_per_time": float(
                analysis["po4_lateral_diffusion_angstrom2_per_time"]
            ),
            "po4_lateral_diffusion_nm2_per_time": float(analysis["po4_lateral_diffusion_nm2_per_time"]),
            "po4_lateral_diffusion_nm2_per_ns": float(analysis["po4_lateral_diffusion_nm2_per_ns"]),
            "po4_lateral_diffusion_um2_per_s": float(analysis["po4_lateral_diffusion_um2_per_s"]),
            "viscosity_proxy_s_per_um2": float(analysis["viscosity_proxy_s_per_um2"]),
            "po4_fit_r2": float(analysis["fit_po4"]["r2"]),
            "n_frames_used": int(analysis["n_frames_used"]),
            "n_po4": int(analysis["n_po4"]),
            "stage7_file": str(result["stage7_file"]),
        }
        if "target_abs_error_um2_per_s" in analysis:
            row["target_abs_error_um2_per_s"] = float(analysis["target_abs_error_um2_per_s"])
        task_rows.append(row)
        grouped.setdefault((row["lj_alpha"], row["slater_alpha"]), []).append(row)

    task_csv = assembled_dir / "task_results.csv"
    with task_csv.open("w", encoding="utf-8", newline="") as fh:
        fieldnames = [
            "analysis_task_id",
            "task_id",
            "code",
            "lj_alpha",
            "slater_alpha",
            "replicate",
            "po4_lateral_diffusion_angstrom2_per_time",
            "po4_lateral_diffusion_nm2_per_time",
            "po4_lateral_diffusion_nm2_per_ns",
            "po4_lateral_diffusion_um2_per_s",
            "viscosity_proxy_s_per_um2",
            "po4_fit_r2",
            "n_frames_used",
            "n_po4",
            "stage7_file",
            "target_abs_error_um2_per_s",
        ]
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in sorted(task_rows, key=lambda item: (item["analysis_task_id"], item["task_id"], item["code"])):
            writer.writerow(row)

    expected_replicates = int(manifest["settings"]["replicates"])
    summary_rows: List[Dict[str, Any]] = []
    for (lj_alpha, slater_alpha), rows in sorted(grouped.items()):
        diffusion_time = np.asarray(
            [row["po4_lateral_diffusion_nm2_per_time"] for row in rows],
            dtype=np.float64,
        )
        diffusion_ns = np.asarray(
            [row["po4_lateral_diffusion_nm2_per_ns"] for row in rows],
            dtype=np.float64,
        )
        diffusion_um = np.asarray(
            [row["po4_lateral_diffusion_um2_per_s"] for row in rows],
            dtype=np.float64,
        )
        viscosity_proxy = np.asarray(
            [row["viscosity_proxy_s_per_um2"] for row in rows],
            dtype=np.float64,
        )
        po4_r2 = np.asarray([row["po4_fit_r2"] for row in rows], dtype=np.float64)
        summary = {
            "lj_alpha": float(lj_alpha),
            "slater_alpha": float(slater_alpha),
            "n_replicates_expected": expected_replicates,
            "n_replicates_completed": len(rows),
            "po4_diffusion_nm2_per_time_mean": float(np.mean(diffusion_time)),
            "po4_diffusion_nm2_per_time_std": float(np.std(diffusion_time, ddof=0)),
            "po4_diffusion_nm2_per_ns_mean": float(np.mean(diffusion_ns)),
            "po4_diffusion_nm2_per_ns_std": float(np.std(diffusion_ns, ddof=0)),
            "po4_diffusion_um2_per_s_mean": float(np.mean(diffusion_um)),
            "po4_diffusion_um2_per_s_std": float(np.std(diffusion_um, ddof=0)),
            "viscosity_proxy_s_per_um2_mean": float(np.mean(viscosity_proxy)),
            "viscosity_proxy_s_per_um2_std": float(np.std(viscosity_proxy, ddof=0)),
            "po4_fit_r2_min": float(np.min(po4_r2)),
        }
        if target_diffusion_um2_s is not None:
            summary["target_diffusion_um2_per_s"] = float(target_diffusion_um2_s)
            summary["target_abs_error_um2_per_s_mean"] = float(
                abs(summary["po4_diffusion_um2_per_s_mean"] - float(target_diffusion_um2_s))
            )
        summary_rows.append(summary)

    summary_csv = assembled_dir / "condition_summary.csv"
    with summary_csv.open("w", encoding="utf-8", newline="") as fh:
        fieldnames = [
            "lj_alpha",
            "slater_alpha",
            "n_replicates_expected",
            "n_replicates_completed",
            "po4_diffusion_nm2_per_time_mean",
            "po4_diffusion_nm2_per_time_std",
            "po4_diffusion_nm2_per_ns_mean",
            "po4_diffusion_nm2_per_ns_std",
            "po4_diffusion_um2_per_s_mean",
            "po4_diffusion_um2_per_s_std",
            "viscosity_proxy_s_per_um2_mean",
            "viscosity_proxy_s_per_um2_std",
            "po4_fit_r2_min",
            "target_diffusion_um2_per_s",
            "target_abs_error_um2_per_s_mean",
        ]
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in summary_rows:
            writer.writerow(row)

    failed_results: List[Dict[str, Any]] = []
    if results_dir.exists():
        for path in sorted(results_dir.glob("*.json")):
            payload = _load_json(path)
            if not payload.get("success"):
                task = payload.get("task", {})
                failed_results.append(
                    {
                        "analysis_task_id": int(payload.get("analysis_task_id", -1)),
                        "task_id": int(task.get("task_id", -1)),
                        "code": str(task.get("code", path.stem)),
                        "error": str(payload.get("error", "")),
                        "stage7_file": str(payload.get("stage7_file", "")),
                    }
                )

    failed_csv = assembled_dir / "failed_tasks.csv"
    with failed_csv.open("w", encoding="utf-8", newline="") as fh:
        fieldnames = ["analysis_task_id", "task_id", "code", "error", "stage7_file"]
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in failed_results:
            writer.writerow(row)

    finite_rows = [row for row in summary_rows if math.isfinite(row["po4_diffusion_um2_per_s_mean"])]
    complete_rows = [
        row for row in finite_rows if int(row["n_replicates_completed"]) == int(row["n_replicates_expected"])
    ]
    recommendation_pool = complete_rows or finite_rows
    recommendation: Dict[str, Any]
    if not recommendation_pool:
        recommendation = {
            "created_at_utc": _now_utc(),
            "base_dir": str(base_dir),
            "selection_basis": "no_finite_conditions",
            "recommended_condition": None,
        }
    else:
        if target_diffusion_um2_s is not None:
            best = min(
                recommendation_pool,
                key=lambda row: (
                    row["target_abs_error_um2_per_s_mean"],
                    row["po4_diffusion_um2_per_s_std"],
                    -row["n_replicates_completed"],
                    row["lj_alpha"],
                    row["slater_alpha"],
                ),
            )
            selection_basis = "closest_target_diffusion_um2_per_s"
        else:
            best = max(
                recommendation_pool,
                key=lambda row: (
                    row["po4_diffusion_um2_per_s_mean"],
                    -row["po4_diffusion_um2_per_s_std"],
                    row["n_replicates_completed"],
                    -row["lj_alpha"],
                    -row["slater_alpha"],
                ),
            )
            selection_basis = "highest_mean_diffusion_um2_per_s"
        recommendation = {
            "created_at_utc": _now_utc(),
            "base_dir": str(base_dir),
            "selection_basis": selection_basis,
            "target_diffusion_um2_per_s": target_diffusion_um2_s,
            "candidate_scope": "complete_replicates" if complete_rows else "finite_results_only",
            "recommended_condition": best,
        }

    recommendation_path = assembled_dir / "recommendation_summary.json"
    _write_json(recommendation_path, recommendation)

    payload = {
        "created_at_utc": _now_utc(),
        "base_dir": str(base_dir),
        "n_analysis_tasks_total": int(len(analysis_manifest["tasks"])),
        "n_analysis_tasks_completed_successfully": len(task_rows),
        "n_analysis_tasks_failed": len(failed_results),
        "task_results_csv": str(task_csv),
        "condition_summary_csv": str(summary_csv),
        "failed_tasks_csv": str(failed_csv),
        "recommendation_summary_json": str(recommendation_path),
    }
    _write_json(assembled_dir / "summary.json", payload)
    print(f"Assembled analysis results written under {assembled_dir}")
    return 0


def cmd_assemble_analysis(args: argparse.Namespace) -> int:
    base_dir = Path(args.base_dir).expanduser().resolve()
    return assemble_analysis(base_dir)


def _array_script_content(round_manifest_path: Path, base_dir: Path, n_tasks: int) -> str:
    output_path = _slurm_dir(base_dir) / "sweep-%A_%a.out"
    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name=hybrid-soften-{base_dir.name}",
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
        f"#SBATCH --job-name=hybrid-soften-collect-{base_dir.name}",
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


def _analysis_array_script_content(round_manifest_path: Path, base_dir: Path, n_tasks: int) -> str:
    output_path = _analysis_slurm_dir(base_dir) / "analyze-%A_%a.out"
    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name=hybrid-analyze-{base_dir.name}",
        f"#SBATCH --output={output_path}",
        f"#SBATCH --time={os.environ.get('HYBRID_SWEEP_ANALYSIS_WALLTIME', '08:00:00')}",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        f"#SBATCH --cpus-per-task={int(os.environ.get('HYBRID_SWEEP_ANALYSIS_CPUS_PER_TASK', '1'))}",
        f"#SBATCH --array=0-{n_tasks - 1}",
    ]
    lines.extend(_optional_sbatch_directives("ANALYSIS"))
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
            'cmd=("$PYTHON_BIN" -u "$SCRIPT_DIR/workflow.py" run-analysis-task --round-manifest "$ROUND_MANIFEST" --task-id "$SLURM_ARRAY_TASK_ID")',
            'if [ "${HYBRID_SWEEP_OVERWRITE:-0}" = "1" ]; then',
            '  cmd+=(--overwrite)',
            "fi",
            '"${cmd[@]}"',
        ]
    )
    return "\n".join(lines) + "\n"


def _analysis_collector_script_content(base_dir: Path) -> str:
    output_path = _analysis_slurm_dir(base_dir) / "collect-analysis-%j.out"
    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name=hybrid-analyze-collect-{base_dir.name}",
        f"#SBATCH --output={output_path}",
        f"#SBATCH --time={os.environ.get('HYBRID_SWEEP_ANALYSIS_COLLECT_WALLTIME', '01:00:00')}",
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
            'export PYTHONUNBUFFERED=1',
            '"$PYTHON_BIN" -u "$SCRIPT_DIR/workflow.py" assemble-analysis --base-dir "$BASE_DIR"',
        ]
    )
    return "\n".join(lines) + "\n"


def cmd_submit_slurm(args: argparse.Namespace) -> int:
    base_dir = Path(args.base_dir).expanduser().resolve()
    manifest = _load_manifest(base_dir)
    n_tasks = len(manifest["tasks"])
    if n_tasks <= 0:
        raise RuntimeError("Sweep manifest contains no tasks")
    round_manifest = {
        "schema": SCHEMA_SLURM_ROUND,
        "created_at_utc": _now_utc(),
        "base_dir": str(base_dir),
        "manifest_path": str(_manifest_path(base_dir)),
        "n_tasks": n_tasks,
    }
    round_manifest_path = _slurm_dir(base_dir) / "round_manifest.json"
    _write_json(round_manifest_path, round_manifest)

    array_script = _slurm_dir(base_dir) / "run_array.sbatch"
    collector_script = _slurm_dir(base_dir) / "collect_results.sbatch"
    _write_text(array_script, _array_script_content(round_manifest_path, base_dir, n_tasks), executable=True)
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
            "run_job_id": train_job,
            "collect_job_id": collect_job,
        },
    )
    print(f"Submitted sweep array job: {train_job}")
    print(f"Submitted collector job: {collect_job}")
    return 0


def cmd_submit_analysis_slurm(args: argparse.Namespace) -> int:
    base_dir = Path(args.base_dir).expanduser().resolve()
    analysis_manifest = _load_analysis_manifest(base_dir)
    n_tasks = len(analysis_manifest["tasks"])
    if n_tasks <= 0:
        raise RuntimeError("Analysis manifest contains no tasks")
    round_manifest = {
        "schema": SCHEMA_ANALYSIS_SLURM,
        "created_at_utc": _now_utc(),
        "base_dir": str(base_dir),
        "analysis_manifest_path": str(_analysis_manifest_path(base_dir)),
        "n_tasks": n_tasks,
    }
    round_manifest_path = _analysis_slurm_dir(base_dir) / "round_manifest.json"
    _write_json(round_manifest_path, round_manifest)

    array_script = _analysis_slurm_dir(base_dir) / "analyze_array.sbatch"
    collector_script = _analysis_slurm_dir(base_dir) / "collect_analysis.sbatch"
    _write_text(
        array_script,
        _analysis_array_script_content(round_manifest_path, base_dir, n_tasks),
        executable=True,
    )
    _write_text(
        collector_script,
        _analysis_collector_script_content(base_dir),
        executable=True,
    )

    print(f"Analysis round manifest: {round_manifest_path}")
    print(f"Analysis array script: {array_script}")
    print(f"Analysis collector script: {collector_script}")

    if args.no_submit:
        print("Staged analysis Slurm scripts only; no jobs submitted.")
        return 0

    if not _shutil_which("sbatch"):
        raise RuntimeError("sbatch not available in PATH")

    analyze_job = _run_sbatch(["sbatch", "--parsable", str(array_script)])
    collect_job = _run_sbatch(["sbatch", "--parsable", f"--dependency=afterok:{analyze_job}", str(collector_script)])
    _write_json(
        _analysis_slurm_dir(base_dir) / "submission.json",
        {
            "created_at_utc": _now_utc(),
            "analyze_job_id": analyze_job,
            "collect_job_id": collect_job,
        },
    )
    print(f"Submitted analysis array job: {analyze_job}")
    print(f"Submitted analysis collector job: {collect_job}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Hybrid interface softening sweep workflow")
    sub = parser.add_subparsers(dest="command", required=True)

    p_init = sub.add_parser("init-run", help="Initialize a new bilayer-only softening sweep")
    p_init.add_argument("--base-dir", default=str(WORKFLOW_DIR / "runs" / "default"))
    p_init.add_argument("--pdb-id", default=DEFAULT_PDB_ID)
    p_init.add_argument("--lj-alphas", type=_float_list_arg, default=DEFAULT_LJ_ALPHAS)
    p_init.add_argument("--slater-alphas", type=_float_list_arg, default=DEFAULT_SLATER_ALPHAS)
    p_init.add_argument("--replicates", type=int, default=DEFAULT_REPLICATES)
    p_init.add_argument("--seed", type=int, default=DEFAULT_SEED)
    p_init.add_argument("--integration-ps-per-step", type=float, default=DEFAULT_INTEGRATION_PS_PER_STEP)
    p_init.add_argument("--target-diffusion-um2-s", type=float, default=None)
    p_init.set_defaults(func=cmd_init_run)

    p_local = sub.add_parser("run-local", help="Run sweep tasks locally")
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

    p_init_analysis = sub.add_parser("init-analysis", help="Discover completed stage-7 files for analysis")
    p_init_analysis.add_argument("--base-dir", default=str(WORKFLOW_DIR / "runs" / "default"))
    p_init_analysis.set_defaults(func=cmd_init_analysis)

    p_analysis_local = sub.add_parser("run-analysis-local", help="Run completed analysis tasks locally")
    p_analysis_local.add_argument("--base-dir", default=str(WORKFLOW_DIR / "runs" / "default"))
    p_analysis_local.add_argument("--max-tasks", type=int, default=None)
    p_analysis_local.add_argument("--start-task", type=int, default=0)
    p_analysis_local.add_argument("--overwrite", action="store_true")
    p_analysis_local.add_argument("--no-assemble", action="store_true")
    p_analysis_local.set_defaults(func=cmd_run_analysis_local)

    p_analysis_task = sub.add_parser("run-analysis-task", help="Run a single analysis Slurm array task")
    p_analysis_task.add_argument("--round-manifest", required=True)
    p_analysis_task.add_argument("--task-id", required=True, type=int)
    p_analysis_task.add_argument("--overwrite", action="store_true")
    p_analysis_task.set_defaults(func=cmd_run_analysis_task)

    p_assemble_analysis = sub.add_parser("assemble-analysis", help="Aggregate completed analysis results")
    p_assemble_analysis.add_argument("--base-dir", default=str(WORKFLOW_DIR / "runs" / "default"))
    p_assemble_analysis.set_defaults(func=cmd_assemble_analysis)

    p_submit_analysis = sub.add_parser("submit-analysis-slurm", help="Stage and optionally submit analysis Slurm scripts")
    p_submit_analysis.add_argument("--base-dir", default=str(WORKFLOW_DIR / "runs" / "default"))
    p_submit_analysis.add_argument("--no-submit", action="store_true")
    p_submit_analysis.set_defaults(func=cmd_submit_analysis_slurm)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
