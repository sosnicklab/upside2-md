#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shutil
import subprocess as sp
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence

import h5py
import numpy as np


WORKFLOW_DIR = Path(os.path.abspath(__file__)).parent
REPO_ROOT = WORKFLOW_DIR.parent
PY_DIR = REPO_ROOT / "py"
if str(PY_DIR) not in sys.path:
    sys.path.insert(0, str(PY_DIR))

import run_upside as ru


PDB_TO_INITIAL_STRUCTURE = PY_DIR / "PDB_to_initial_structure.py"
UPSIDE_EXECUTABLE = REPO_ROOT / "obj" / "upside"
HYBRID_RUN_SCRIPT = REPO_ROOT / "example" / "16.MARTINI" / "run_sim_1rkl.sh"
REFERENCE_PDB_DIR = REPO_ROOT / "example" / "08.MembraneSimulation" / "pdb"
HYBRID_PDB_DIR = REPO_ROOT / "example" / "16.MARTINI" / "pdb"

SCHEMA_MANIFEST = "hybrid_interface_rmsf_sweep_manifest_v1"
SCHEMA_SLURM_ROUND = "hybrid_interface_rmsf_sweep_slurm_round_v1"
SCHEMA_TASK_RESULT = "hybrid_interface_rmsf_sweep_task_result_v1"
SCHEMA_ANALYSIS_MANIFEST = "hybrid_interface_rmsf_sweep_analysis_manifest_v2"
SCHEMA_ANALYSIS_SLURM = "hybrid_interface_rmsf_sweep_analysis_slurm_v2"
SCHEMA_ANALYSIS_RESULT = "hybrid_interface_rmsf_sweep_analysis_result_v2"

REFERENCE_METHOD = "example_08_fixed_curvature"
HYBRID_METHOD = "example_16_hybrid"

DEFAULT_PDB_ID = "1rkl"
DEFAULT_INTERFACE_SCALES = [
    0.825,
    0.85,
    0.875,
    0.90,
    0.925,
    0.95,
    0.975,
    1.00,
    1.025,
    1.05,
    1.075,
    1.10,
    1.125,
    1.15,
    1.175,
    1.20,
]
DEFAULT_HYBRID_REPLICATES = 1
DEFAULT_REFERENCE_REPLICATES = 4
DEFAULT_SEED = 20260414
DEFAULT_BURN_IN_FRACTION = 0.20
DEFAULT_TRENDLINE_SAMPLES = 201
DEFAULT_STABILITY_MEAN_RMSF_RATIO_MAX = 3.0
DEFAULT_STABILITY_MAX_RMSF_RATIO_MAX = 3.0
DEFAULT_STABILITY_CA_RG_RATIO_MAX = 1.75
DEFAULT_STABILITY_CA_SPAN_RATIO_MAX = 1.75
DEFAULT_EMBEDDED_OCCUPANCY_MIN = 0.50

DEFAULT_REFERENCE_SETTINGS = {
    "REFERENCE_TEMPERATURE": "0.80",
    "REFERENCE_DURATION": "200001",
    "REFERENCE_FRAME_INTERVAL": "100",
    "REFERENCE_MEMBRANE_THICKNESS": "24.8",
    "REFERENCE_USE_CURVATURE": "1",
    "REFERENCE_CURVATURE_RADIUS": "120.0",
    "REFERENCE_CURVATURE_SIGN": "1",
}

DEFAULT_HYBRID_RUNTIME_SETTINGS = {
    "EQ_62_NSTEPS": "1000",
    "EQ_63_NSTEPS": "1000",
    "EQ_64_NSTEPS": "1000",
    "EQ_65_NSTEPS": "1000",
    "EQ_66_NSTEPS": "1000",
    "PROD_70_NSTEPS": "50000",
    "EQ_FRAME_STEPS": "250",
    "PROD_FRAME_STEPS": "100",
}

DEFAULT_REFERENCE_SETTING_KEYS = [
    "REFERENCE_TEMPERATURE",
    "REFERENCE_DURATION",
    "REFERENCE_FRAME_INTERVAL",
    "REFERENCE_MEMBRANE_THICKNESS",
    "REFERENCE_USE_CURVATURE",
    "REFERENCE_CURVATURE_RADIUS",
    "REFERENCE_CURVATURE_SIGN",
]

DEFAULT_HYBRID_PASSTHROUGH_ENV_KEYS = [
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
    hybrid_replicates: int
    reference_replicates: int
    seed: int
    burn_in_fraction: float
    trendline_samples: int
    embedded_occupancy_min: float
    stability_mean_rmsf_ratio_max: float
    stability_max_rmsf_ratio_max: float
    stability_ca_rg_ratio_max: float
    stability_ca_span_ratio_max: float
    reference_settings: Dict[str, str]
    hybrid_passthrough_env: Dict[str, str]


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


def _truthy(text: str | None) -> bool:
    if text is None:
        return False
    return text.strip().lower() not in {"", "0", "false", "no", "off"}


def _capture_reference_settings() -> Dict[str, str]:
    extra = _floatless_csv(os.environ.get("HYBRID_SWEEP_EXTRA_REFERENCE_KEYS", ""))
    keys = _normalize_env_key_list(DEFAULT_REFERENCE_SETTING_KEYS + extra)
    out = dict(DEFAULT_REFERENCE_SETTINGS)
    for key in keys:
        value = os.environ.get(key)
        if value is None or value == "":
            continue
        out[key] = value
    return out


def _capture_hybrid_passthrough_env() -> Dict[str, str]:
    extra = _floatless_csv(os.environ.get("HYBRID_SWEEP_EXTRA_ENV_KEYS", ""))
    keys = _normalize_env_key_list(DEFAULT_HYBRID_PASSTHROUGH_ENV_KEYS + extra)
    out: Dict[str, str] = dict(DEFAULT_HYBRID_RUNTIME_SETTINGS)
    for key in keys:
        value = os.environ.get(key)
        if value is None or value == "":
            continue
        out[key] = value
    return out


def _build_config(args: argparse.Namespace, base_dir: Path) -> Config:
    interface_scales = sorted(dict.fromkeys(float(x) for x in args.interface_scales))
    for value in interface_scales:
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"interface scale must be finite and > 0, got {value}")
    hybrid_replicates = int(args.hybrid_replicates)
    reference_replicates = int(args.reference_replicates)
    if hybrid_replicates < 1:
        raise ValueError(f"hybrid replicates must be >= 1, got {hybrid_replicates}")
    if reference_replicates < 1:
        raise ValueError(f"reference replicates must be >= 1, got {reference_replicates}")
    burn_in_fraction = float(args.burn_in_fraction)
    if not math.isfinite(burn_in_fraction) or not (0.0 <= burn_in_fraction < 1.0):
        raise ValueError(f"burn-in fraction must be finite in [0, 1), got {burn_in_fraction}")
    trendline_samples = int(args.trendline_samples)
    if trendline_samples < 11:
        raise ValueError(f"trendline samples must be >= 11, got {trendline_samples}")
    embedded_occupancy_min = float(
        os.environ.get(
            "HYBRID_SWEEP_EMBEDDED_OCCUPANCY_MIN",
            str(DEFAULT_EMBEDDED_OCCUPANCY_MIN),
        )
    )
    if not math.isfinite(embedded_occupancy_min) or not (0.0 < embedded_occupancy_min <= 1.0):
        raise ValueError(
            f"embedded occupancy min must be finite in (0, 1], got {embedded_occupancy_min}"
        )
    stability_mean_rmsf_ratio_max = float(
        os.environ.get(
            "HYBRID_SWEEP_STABILITY_MEAN_RMSF_RATIO_MAX",
            str(DEFAULT_STABILITY_MEAN_RMSF_RATIO_MAX),
        )
    )
    stability_max_rmsf_ratio_max = float(
        os.environ.get(
            "HYBRID_SWEEP_STABILITY_MAX_RMSF_RATIO_MAX",
            str(DEFAULT_STABILITY_MAX_RMSF_RATIO_MAX),
        )
    )
    stability_ca_rg_ratio_max = float(
        os.environ.get(
            "HYBRID_SWEEP_STABILITY_CA_RG_RATIO_MAX",
            str(DEFAULT_STABILITY_CA_RG_RATIO_MAX),
        )
    )
    stability_ca_span_ratio_max = float(
        os.environ.get(
            "HYBRID_SWEEP_STABILITY_CA_SPAN_RATIO_MAX",
            str(DEFAULT_STABILITY_CA_SPAN_RATIO_MAX),
        )
    )
    for name, value in (
        ("stability_mean_rmsf_ratio_max", stability_mean_rmsf_ratio_max),
        ("stability_max_rmsf_ratio_max", stability_max_rmsf_ratio_max),
        ("stability_ca_rg_ratio_max", stability_ca_rg_ratio_max),
        ("stability_ca_span_ratio_max", stability_ca_span_ratio_max),
    ):
        if not math.isfinite(value) or value <= 0.0:
            raise ValueError(f"{name} must be finite and > 0, got {value}")
    return Config(
        base_dir=base_dir,
        pdb_id=str(args.pdb_id),
        interface_scales=interface_scales,
        hybrid_replicates=hybrid_replicates,
        reference_replicates=reference_replicates,
        seed=int(args.seed),
        burn_in_fraction=burn_in_fraction,
        trendline_samples=trendline_samples,
        embedded_occupancy_min=embedded_occupancy_min,
        stability_mean_rmsf_ratio_max=stability_mean_rmsf_ratio_max,
        stability_max_rmsf_ratio_max=stability_max_rmsf_ratio_max,
        stability_ca_rg_ratio_max=stability_ca_rg_ratio_max,
        stability_ca_span_ratio_max=stability_ca_span_ratio_max,
        reference_settings=_capture_reference_settings(),
        hybrid_passthrough_env=_capture_hybrid_passthrough_env(),
    )


def _build_manifest(config: Config) -> Dict[str, Any]:
    tasks: List[Dict[str, Any]] = []
    task_id = 0
    for replicate in range(config.reference_replicates):
        code = f"ref_r{replicate + 1:02d}"
        tasks.append(
            {
                "task_id": task_id,
                "kind": "reference",
                "method": REFERENCE_METHOD,
                "code": code,
                "replicate": int(replicate + 1),
                "seed": _task_seed(config.seed, task_id),
            }
        )
        task_id += 1
    for interface_scale in config.interface_scales:
        for replicate in range(config.hybrid_replicates):
            code = f"scale{_format_float_tag(interface_scale)}_r{replicate + 1:02d}"
            tasks.append(
                {
                    "task_id": task_id,
                    "kind": "hybrid",
                    "method": HYBRID_METHOD,
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
            "hybrid_replicates": config.hybrid_replicates,
            "reference_replicates": config.reference_replicates,
            "seed": config.seed,
            "reference_settings": config.reference_settings,
            "hybrid_passthrough_env": config.hybrid_passthrough_env,
            "reference_pdb_dir": str(REFERENCE_PDB_DIR),
            "hybrid_pdb_dir": str(HYBRID_PDB_DIR),
            "pdb_to_initial_structure": str(PDB_TO_INITIAL_STRUCTURE),
            "upside_executable": str(UPSIDE_EXECUTABLE),
            "hybrid_run_script": str(HYBRID_RUN_SCRIPT),
            "analysis_settings": {
                "burn_in_fraction": config.burn_in_fraction,
                "trendline_samples": config.trendline_samples,
                "embedded_occupancy_min": config.embedded_occupancy_min,
                "rmsf_atom_role": "CA",
                "alignment_iterations": 3,
                "stability_mean_rmsf_ratio_max": config.stability_mean_rmsf_ratio_max,
                "stability_max_rmsf_ratio_max": config.stability_max_rmsf_ratio_max,
                "stability_ca_rg_ratio_max": config.stability_ca_rg_ratio_max,
                "stability_ca_span_ratio_max": config.stability_ca_span_ratio_max,
            },
        },
        "tasks": tasks,
    }


def _reference_pdb_path(pdb_id: str) -> Path:
    return REFERENCE_PDB_DIR / f"{pdb_id}.pdb"


def _hybrid_pdb_path(pdb_id: str) -> Path:
    return HYBRID_PDB_DIR / f"{pdb_id}.pdb"


def _require_runtime_files(pdb_id: str) -> None:
    if not PDB_TO_INITIAL_STRUCTURE.exists():
        raise RuntimeError(f"PDB-to-initial-structure script not found: {PDB_TO_INITIAL_STRUCTURE}")
    if not UPSIDE_EXECUTABLE.exists():
        raise RuntimeError(f"Upside executable not found: {UPSIDE_EXECUTABLE}")
    if not HYBRID_RUN_SCRIPT.exists():
        raise RuntimeError(f"Hybrid run script not found: {HYBRID_RUN_SCRIPT}")
    reference_pdb = _reference_pdb_path(pdb_id)
    if not reference_pdb.exists():
        raise RuntimeError(f"Reference PDB not found: {reference_pdb}")
    hybrid_pdb = _hybrid_pdb_path(pdb_id)
    if not hybrid_pdb.exists():
        raise RuntimeError(f"Hybrid PDB not found: {hybrid_pdb}")


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
    if (REPO_ROOT / ".venv" / "bin").exists():
        _prepend_env_path(env, "PATH", str(REPO_ROOT / ".venv" / "bin"))
    _prepend_env_path(env, "PATH", str(REPO_ROOT / "obj"))
    _prepend_env_path(env, "PYTHONPATH", str(REPO_ROOT / "py"))
    return env


def _call_logged(cmd: Sequence[str], log_path: Path, env: Dict[str, str], cwd: Path | None = None) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as fh:
        proc = sp.run(
            list(cmd),
            cwd=str(cwd if cwd is not None else REPO_ROOT),
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


def _write_log_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _reference_config_kwargs(pdb_id: str, input_dir: Path, settings: Dict[str, str]) -> Dict[str, Any]:
    param_dir_common = REPO_ROOT / "parameters" / "common"
    param_dir_ff = REPO_ROOT / "parameters" / "ff_2.1"
    kwargs: Dict[str, Any] = {
        "rama_library": str(param_dir_common / "rama.dat"),
        "rama_sheet_mix_energy": str(param_dir_ff / "sheet"),
        "reference_state_rama": str(param_dir_common / "rama_reference.pkl"),
        "hbond_energy": str(param_dir_ff / "hbond.h5"),
        "rotamer_placement": str(param_dir_ff / "sidechain.h5"),
        "dynamic_rotamer_1body": True,
        "rotamer_interaction": str(param_dir_ff / "sidechain.h5"),
        "environment_potential": str(param_dir_ff / "environment.h5"),
        "bb_environment_potential": str(param_dir_ff / "bb_env.dat"),
        "membrane_potential": str(param_dir_ff / "membrane.h5"),
        "membrane_thickness": float(settings["REFERENCE_MEMBRANE_THICKNESS"]),
        "chain_break_from_file": str(input_dir / f"{pdb_id}.chain_breaks"),
        "initial_structure": str(input_dir / f"{pdb_id}.initial.npy"),
    }
    if _truthy(settings.get("REFERENCE_USE_CURVATURE")):
        kwargs["use_curvature"] = True
        kwargs["curvature_radius"] = float(settings["REFERENCE_CURVATURE_RADIUS"])
        kwargs["curvature_sign"] = int(settings["REFERENCE_CURVATURE_SIGN"])
    return kwargs


def _prepare_reference_inputs(
    *,
    pdb_id: str,
    input_dir: Path,
    settings: Dict[str, str],
    env: Dict[str, str],
    prep_log: Path,
    config_log: Path,
) -> Path:
    reference_pdb = _reference_pdb_path(pdb_id)
    input_dir.mkdir(parents=True, exist_ok=True)
    prefix = input_dir / pdb_id
    cmd = [
        sys.executable,
        str(PDB_TO_INITIAL_STRUCTURE),
        str(reference_pdb),
        str(prefix),
        "--record-chain-breaks",
        "--disable-recentering",
    ]
    _call_logged(cmd, prep_log, env)

    fasta = input_dir / f"{pdb_id}.fasta"
    config_base = input_dir / f"{pdb_id}.reference.base.up"
    kwargs = _reference_config_kwargs(pdb_id, input_dir, settings)
    output = ru.upside_config(str(fasta), str(config_base), **kwargs)
    _write_log_text(config_log, output)
    return config_base


def _run_reference_md(
    *,
    output_file: Path,
    settings: Dict[str, str],
    seed: int,
    env: Dict[str, str],
    log_path: Path,
) -> Dict[str, Any]:
    cmd = [
        str(UPSIDE_EXECUTABLE),
        str(output_file),
        "--duration",
        str(settings["REFERENCE_DURATION"]),
        "--frame-interval",
        str(settings["REFERENCE_FRAME_INTERVAL"]),
        "--temperature",
        str(settings["REFERENCE_TEMPERATURE"]),
        "--seed",
        str(int(seed)),
        "--disable-recentering",
        "--integrator",
        "mv",
    ]
    _call_logged(cmd, log_path, env)
    return {
        "duration": float(settings["REFERENCE_DURATION"]),
        "frame_interval": float(settings["REFERENCE_FRAME_INTERVAL"]),
        "temperature": float(settings["REFERENCE_TEMPERATURE"]),
        "integrator": "mv",
        "disable_recentering": True,
    }


def _run_reference_task(
    *,
    base_dir: Path,
    manifest: Dict[str, Any],
    task: Dict[str, Any],
) -> Dict[str, Any]:
    settings = manifest["settings"]
    reference_settings = dict(settings["reference_settings"])
    task_dir = _task_dir(base_dir, task["code"])
    input_dir = task_dir / "inputs"
    run_dir = task_dir / "run"
    log_dir = task_dir / "logs"
    run_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    env = _base_runtime_env()
    env["SEED"] = str(int(task["seed"]))
    env["OMP_NUM_THREADS"] = os.environ.get("OMP_NUM_THREADS", "1")

    prep_log = log_dir / "prepare_reference_inputs.log"
    config_log = log_dir / "configure_reference.log"
    run_log = log_dir / "run_reference.log"

    config_base = _prepare_reference_inputs(
        pdb_id=str(settings["pdb_id"]),
        input_dir=input_dir,
        settings=reference_settings,
        env=env,
        prep_log=prep_log,
        config_log=config_log,
    )

    output_file = run_dir / f"{settings['pdb_id']}.reference.r{int(task['replicate']):02d}.up"
    shutil.copy2(config_base, output_file)
    runtime = _run_reference_md(
        output_file=output_file,
        settings=reference_settings,
        seed=int(task["seed"]),
        env=env,
        log_path=run_log,
    )
    if not output_file.exists():
        raise RuntimeError(f"Reference run did not produce output file: {output_file}")

    return {
        "schema": SCHEMA_TASK_RESULT,
        "created_at_utc": _now_utc(),
        "base_dir": str(base_dir),
        "task": task,
        "success": True,
        "kind": "reference",
        "method": REFERENCE_METHOD,
        "run_dir": str(run_dir),
        "log_dir": str(log_dir),
        "analysis_target_file": str(output_file),
        "log_files": {
            "prepare_reference_inputs": str(prep_log),
            "configure_reference": str(config_log),
            "run_reference": str(run_log),
        },
        "reference_runtime": runtime,
    }


def _run_hybrid_task(
    *,
    base_dir: Path,
    manifest: Dict[str, Any],
    task: Dict[str, Any],
) -> Dict[str, Any]:
    settings = manifest["settings"]
    passthrough_env = dict(settings["hybrid_passthrough_env"])
    task_dir = _task_dir(base_dir, task["code"])
    run_dir = task_dir / "run"
    log_dir = task_dir / "logs"
    run_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    env = _base_runtime_env()
    env.update({str(k): str(v) for k, v in passthrough_env.items()})
    env["SEED"] = str(int(task["seed"]))
    env.setdefault("PREP_SEED", env["SEED"])
    env["RUN_DIR"] = str(run_dir)
    env["PROTEIN_ENV_INTERFACE_SCALE"] = f"{float(task['interface_scale']):.10g}"
    env.setdefault("PROTEIN_AA_PDB", str(_hybrid_pdb_path(str(settings["pdb_id"]))))
    env["HYBRID_SWEEP_PROJECT_ROOT"] = str(REPO_ROOT)

    run_log = log_dir / "run_hybrid_workflow.log"
    cmd = ["bash", str(HYBRID_RUN_SCRIPT), f"PDB_ID={settings['pdb_id']}"]
    _call_logged(cmd, run_log, env, cwd=HYBRID_RUN_SCRIPT.parent)

    stage_70_file = run_dir / "checkpoints" / f"{settings['pdb_id']}.stage_7.0.up"
    if not stage_70_file.exists():
        raise RuntimeError(f"Hybrid workflow did not produce stage-7 checkpoint: {stage_70_file}")

    return {
        "schema": SCHEMA_TASK_RESULT,
        "created_at_utc": _now_utc(),
        "base_dir": str(base_dir),
        "task": task,
        "success": True,
        "kind": "hybrid",
        "method": HYBRID_METHOD,
        "run_dir": str(run_dir),
        "log_dir": str(log_dir),
        "analysis_target_file": str(stage_70_file),
        "log_files": {
            "run_hybrid_workflow": str(run_log),
        },
        "hybrid_runtime": {
            "interface_scale": float(task["interface_scale"]),
        },
    }


def cmd_init_run(args: argparse.Namespace) -> int:
    base_dir = Path(args.base_dir).expanduser().resolve()
    base_dir.mkdir(parents=True, exist_ok=True)
    config = _build_config(args, base_dir)
    _require_runtime_files(config.pdb_id)
    manifest = _build_manifest(config)
    _write_json(_manifest_path(base_dir), manifest)
    print(f"Initialized hybrid interface RMSF sweep: {base_dir}")
    print(f"Reference tasks: {config.reference_replicates}")
    print(f"Hybrid tasks: {len(config.interface_scales) * config.hybrid_replicates}")
    print(f"Interface scales: {manifest['settings']['interface_scales']}")
    return 0


def _run_task(base_dir: Path, manifest: Dict[str, Any], task: Dict[str, Any], overwrite: bool) -> Dict[str, Any]:
    result_path = _result_path(base_dir, task["code"])
    if result_path.exists() and not overwrite:
        existing = _load_json(result_path)
        if existing.get("success"):
            print(f"Skipping completed task {task['code']}")
            return existing

    try:
        if task["kind"] == "reference":
            result = _run_reference_task(base_dir=base_dir, manifest=manifest, task=task)
        elif task["kind"] == "hybrid":
            result = _run_hybrid_task(base_dir=base_dir, manifest=manifest, task=task)
        else:
            raise RuntimeError(f"Unsupported task kind: {task['kind']}")
        _write_json(result_path, result)
        return result
    except Exception as exc:
        result = {
            "schema": SCHEMA_TASK_RESULT,
            "created_at_utc": _now_utc(),
            "base_dir": str(base_dir),
            "task": task,
            "success": False,
            "kind": str(task.get("kind", "")),
            "method": str(task.get("method", "")),
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


def _successful_analysis_results_from_dir(result_dir: Path) -> List[Dict[str, Any]]:
    results: List[Dict[str, Any]] = []
    if not result_dir.exists():
        return results
    for path in sorted(result_dir.glob("*.json")):
        payload = _load_json(path)
        if payload.get("schema") != SCHEMA_ANALYSIS_RESULT:
            continue
        if payload.get("success"):
            results.append(payload)
    return results


def assemble_results(base_dir: Path) -> int:
    manifest = _load_manifest(base_dir)
    results = _successful_results_from_dir(_result_dir(base_dir))
    assembled_dir = _assembled_dir(base_dir)
    assembled_dir.mkdir(parents=True, exist_ok=True)

    task_rows: List[Dict[str, Any]] = []
    hybrid_grouped: Dict[float, List[Dict[str, Any]]] = {}
    n_reference_completed = 0
    for result in results:
        task = result["task"]
        row = {
            "task_id": int(task["task_id"]),
            "code": str(task["code"]),
            "kind": str(task["kind"]),
            "method": str(task["method"]),
            "interface_scale": (
                ""
                if task["kind"] != "hybrid"
                else f"{float(task['interface_scale']):.10g}"
            ),
            "replicate": int(task["replicate"]),
            "seed": int(task["seed"]),
            "run_dir": str(result["run_dir"]),
            "log_dir": str(result["log_dir"]),
            "analysis_target_file": str(result["analysis_target_file"]),
        }
        task_rows.append(row)
        if task["kind"] == "reference":
            n_reference_completed += 1
        else:
            hybrid_grouped.setdefault(float(task["interface_scale"]), []).append(row)

    task_csv = assembled_dir / "task_results.csv"
    with task_csv.open("w", encoding="utf-8", newline="") as fh:
        fieldnames = [
            "task_id",
            "code",
            "kind",
            "method",
            "interface_scale",
            "replicate",
            "seed",
            "run_dir",
            "log_dir",
            "analysis_target_file",
        ]
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in sorted(task_rows, key=lambda item: item["task_id"]):
            writer.writerow(row)

    summary_rows: List[Dict[str, Any]] = []
    expected_hybrid_replicates = int(manifest["settings"]["hybrid_replicates"])
    for interface_scale, rows in sorted(hybrid_grouped.items()):
        summary_rows.append(
            {
                "interface_scale": float(interface_scale),
                "n_replicates_expected": expected_hybrid_replicates,
                "n_replicates_completed": len(rows),
                "all_outputs_present": int(all(Path(row["analysis_target_file"]).exists() for row in rows)),
            }
        )

    summary_csv = assembled_dir / "condition_summary.csv"
    with summary_csv.open("w", encoding="utf-8", newline="") as fh:
        fieldnames = [
            "interface_scale",
            "n_replicates_expected",
            "n_replicates_completed",
            "all_outputs_present",
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
        "n_reference_tasks_completed": int(n_reference_completed),
        "n_hybrid_tasks_completed": int(sum(len(rows) for rows in hybrid_grouped.values())),
        "n_hybrid_conditions_completed": len(summary_rows),
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


def _decode_str_values(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values).reshape(-1)
    out: List[str] = []
    for item in arr:
        if isinstance(item, (bytes, np.bytes_)):
            out.append(item.decode("utf-8", errors="ignore").strip())
        else:
            out.append(str(item).strip())
    return np.asarray(out, dtype=object)


def _unique_selected_strings(values: np.ndarray, indices: np.ndarray) -> List[str]:
    selected = _decode_str_values(np.asarray(values)[np.asarray(indices, dtype=np.int64)])
    return sorted({str(item) for item in selected})


def _load_box_frames_optional(h5: h5py.File, n_frame: int) -> np.ndarray | None:
    if "output/box" in h5:
        box = _normalize_box(h5["output/box"][:])
        if box.shape[0] != n_frame:
            raise ValueError("output/box length does not match output/pos frames")
        return box
    if "input/potential/martini_potential" in h5:
        grp = h5["input/potential/martini_potential"]
        required = ("x_len", "y_len", "z_len")
        if all(key in grp.attrs for key in required):
            box = np.array([grp.attrs["x_len"], grp.attrs["y_len"], grp.attrs["z_len"]], dtype=np.float64)
            return np.repeat(box[None, :], n_frame, axis=0)
    return None


def _minimum_image(delta: np.ndarray, box: np.ndarray) -> np.ndarray:
    out = np.asarray(delta, dtype=np.float64).copy()
    box = np.asarray(box, dtype=np.float64)
    valid = np.isfinite(box) & (box > 0.0)
    if np.any(valid):
        out[..., valid] -= box[valid] * np.round(out[..., valid] / box[valid])
    return out


def _unwrap_xyz(positions: np.ndarray, box: np.ndarray) -> np.ndarray:
    unwrapped = np.zeros_like(positions, dtype=np.float64)
    unwrapped[0] = positions[0]
    for frame in range(1, positions.shape[0]):
        current = positions[frame]
        previous = positions[frame - 1]
        if not (np.all(np.isfinite(current)) and np.all(np.isfinite(previous))):
            unwrapped[frame] = current
            continue
        delta = current - previous
        delta = _minimum_image(delta, box[frame - 1])
        unwrapped[frame] = unwrapped[frame - 1] + delta
    return unwrapped


def _align_points_kabsch(mobile: np.ndarray, target: np.ndarray) -> np.ndarray:
    if not (np.all(np.isfinite(mobile)) and np.all(np.isfinite(target))):
        raise ValueError("Kabsch alignment received non-finite coordinates")
    mobile_centroid = np.mean(mobile, axis=0)
    target_centroid = np.mean(target, axis=0)
    mobile_centered = mobile - mobile_centroid
    target_centered = target - target_centroid
    cov = mobile_centered.T @ target_centered
    if not np.all(np.isfinite(cov)):
        raise ValueError("Kabsch covariance is non-finite")
    try:
        u, _, vt = np.linalg.svd(cov)
    except np.linalg.LinAlgError:
        jitter = np.eye(3, dtype=np.float64) * 1.0e-12
        u, _, vt = np.linalg.svd(cov + jitter)
    correction = np.eye(3, dtype=np.float64)
    if np.linalg.det(u @ vt) < 0.0:
        correction[-1, -1] = -1.0
    rotation = u @ correction @ vt
    return mobile_centered @ rotation + target_centroid


def _iterative_align_frames(frames: np.ndarray, n_iter: int = 3) -> tuple[np.ndarray, np.ndarray]:
    aligned = np.asarray(frames, dtype=np.float64)
    if aligned.ndim != 3 or aligned.shape[0] < 1:
        raise ValueError(f"Unexpected frame array for alignment: {aligned.shape}")
    if not np.all(np.isfinite(aligned[0])):
        raise ValueError("Alignment reference frame is non-finite")
    reference = np.asarray(aligned[0], dtype=np.float64)
    for _ in range(max(1, int(n_iter))):
        next_aligned = np.empty_like(aligned)
        for index, frame in enumerate(aligned):
            next_aligned[index] = _align_points_kabsch(np.asarray(frame, dtype=np.float64), reference)
        aligned = next_aligned
        reference = np.mean(aligned, axis=0)
    return aligned, reference


def _extract_backbone_map(h5: h5py.File) -> Dict[str, Any]:
    inp = h5["input"]
    residue_labels: np.ndarray | None = None
    if "sequence" in inp:
        residue_labels = _decode_str_values(inp["sequence"][:])

    if "hybrid_bb_map" in inp and "aa_atom_index" in inp["hybrid_bb_map"]:
        atom_index = np.asarray(inp["hybrid_bb_map"]["aa_atom_index"][:], dtype=np.int64)
        source = "input/hybrid_bb_map/aa_atom_index"
    elif "potential" in inp and "affine_alignment" in inp["potential"] and "atoms" in inp["potential"]["affine_alignment"]:
        atom_index = np.asarray(inp["potential"]["affine_alignment"]["atoms"][:], dtype=np.int64)
        source = "input/potential/affine_alignment/atoms"
    else:
        raise ValueError("Could not locate a residue-level backbone map for RMSF analysis")

    if atom_index.ndim != 2 or atom_index.shape[0] < 1 or atom_index.shape[1] < 2:
        raise ValueError(f"Unexpected backbone map shape: {atom_index.shape}")
    if np.any(atom_index < 0):
        raise ValueError("Backbone map contains negative atom indices")

    n_residue = int(atom_index.shape[0])
    if residue_labels is None or residue_labels.shape[0] != n_residue:
        residue_labels = np.asarray([str(index + 1) for index in range(n_residue)], dtype=object)

    return {
        "atom_index": atom_index,
        "source": source,
        "residue_labels": residue_labels,
        "ca_column": 1,
    }


def _extract_reference_membrane_geometry(h5: h5py.File) -> Dict[str, Any]:
    potential = h5["input"]["potential"]
    group_name = ""
    if "cb_membrane_potential" in potential:
        group_name = "cb_membrane_potential"
    elif "cb_surf_membrane_potential" in potential:
        group_name = "cb_surf_membrane_potential"
    else:
        raise ValueError("Could not locate a reference membrane potential group for embedded-region selection")

    grp = potential[group_name]
    if "cb_index" not in grp:
        raise ValueError(f"Membrane potential group missing cb_index: input/potential/{group_name}")
    cb_index = np.asarray(grp["cb_index"][:], dtype=np.int64).reshape(-1)
    if cb_index.size < 1:
        raise ValueError(f"Membrane potential group has no residue indices: input/potential/{group_name}")
    if np.any(cb_index < 0):
        raise ValueError(f"Membrane potential group has negative residue indices: input/potential/{group_name}")
    return {
        "group_path": f"input/potential/{group_name}",
        "half_thickness": float(grp.attrs["half_thickness"]),
        "use_curvature": bool(int(grp.attrs.get("use_curvature", 0))),
        "curvature_radius": float(grp.attrs.get("curvature_radius", 0.0)),
        "curvature_sign": float(grp.attrs.get("curvature_sign", 1.0)),
        "cb_index": cb_index,
    }


def _compute_membrane_depth_for_positions(
    residue_xyz: np.ndarray,
    *,
    half_thickness: float,
    use_curvature: bool,
    curvature_radius: float,
    curvature_sign: float,
) -> np.ndarray:
    xyz = np.asarray(residue_xyz, dtype=np.float64)
    if xyz.ndim != 3 or xyz.shape[-1] != 3:
        raise ValueError(f"Unexpected residue coordinate shape for membrane depth: {xyz.shape}")
    if use_curvature:
        center = np.array([0.0, 0.0, -curvature_radius * curvature_sign], dtype=np.float64)
        dist = np.linalg.norm(xyz - center[None, None, :], axis=2)
        return curvature_sign * (dist - curvature_radius)
    return xyz[:, :, 2]


def _derive_reference_embedded_region(
    tasks: Sequence[Dict[str, Any]],
    *,
    burn_in_fraction: float,
    occupancy_min: float,
) -> Dict[str, Any]:
    reference_tasks = [item for item in tasks if item["task"]["kind"] == "reference"]
    if not reference_tasks:
        raise RuntimeError("Need at least one reference task to derive the embedded residue region")

    reference_labels: List[str] | None = None
    reference_numbers: List[int] | None = None
    occupancy_accum: np.ndarray | None = None
    eligible_mask: np.ndarray | None = None
    geometry_path = ""
    half_thickness = None
    use_curvature = None
    curvature_radius = None
    curvature_sign = None
    source_count = 0

    for item in reference_tasks:
        stage_file = Path(item["analysis_target_file"]).expanduser().resolve()
        if not stage_file.exists():
            raise RuntimeError(f"Reference analysis target file not found for embedded-region selection: {stage_file}")
        with h5py.File(stage_file, "r") as h5:
            if "output/pos" not in h5:
                raise ValueError(f"Missing output/pos in {stage_file}")
            pos = _normalize_output_positions(h5["output/pos"][:])
            n_frame_total = int(pos.shape[0])
            if n_frame_total < 5:
                raise ValueError(
                    f"Need at least 5 output frames to derive the embedded region, found {n_frame_total}"
                )
            burn_start = int(math.floor(float(burn_in_fraction) * n_frame_total))
            burn_start = min(max(0, burn_start), n_frame_total - 1)
            pos = pos[burn_start:]
            if pos.shape[0] < 5:
                raise ValueError("Too few frames remain after burn-in to derive the embedded region")

            backbone = _extract_backbone_map(h5)
            residue_labels = [str(x) for x in backbone["residue_labels"]]
            residue_numbers = [int(idx + 1) for idx in range(len(residue_labels))]
            atom_index = np.asarray(backbone["atom_index"], dtype=np.int64)
            ca_index = np.asarray(atom_index[:, int(backbone["ca_column"])], dtype=np.int64)
            if np.max(ca_index) >= pos.shape[1]:
                raise ValueError("Backbone CA selection refers to atoms outside output/pos")

            geometry = _extract_reference_membrane_geometry(h5)
            geometry_path = str(geometry["group_path"])
            half_thickness = float(geometry["half_thickness"])
            use_curvature = bool(geometry["use_curvature"])
            curvature_radius = float(geometry["curvature_radius"])
            curvature_sign = float(geometry["curvature_sign"])

            ca_xyz = np.asarray(pos[:, ca_index, :], dtype=np.float64)
            finite_mask = np.all(np.isfinite(ca_xyz), axis=(1, 2))
            ca_xyz = ca_xyz[finite_mask]
            if ca_xyz.shape[0] < 5:
                raise ValueError(
                    f"Too few finite CA frames remain to derive the embedded region from {stage_file}: {ca_xyz.shape[0]}"
                )
            depth = _compute_membrane_depth_for_positions(
                ca_xyz,
                half_thickness=float(geometry["half_thickness"]),
                use_curvature=bool(geometry["use_curvature"]),
                curvature_radius=float(geometry["curvature_radius"]),
                curvature_sign=float(geometry["curvature_sign"]),
            )
            current_eligible = np.zeros(len(residue_labels), dtype=bool)
            current_eligible[np.asarray(geometry["cb_index"], dtype=np.int64)] = True
            occupancy = np.mean(np.abs(depth) <= float(geometry["half_thickness"]), axis=0)
            occupancy = np.where(current_eligible, occupancy, 0.0)

            if reference_labels is None:
                reference_labels = residue_labels
                reference_numbers = residue_numbers
                occupancy_accum = np.asarray(occupancy, dtype=np.float64)
                eligible_mask = np.asarray(current_eligible, dtype=bool)
            else:
                if residue_labels != reference_labels:
                    raise RuntimeError("Residue labels differ across reference tasks while deriving the embedded region")
                if residue_numbers != reference_numbers:
                    raise RuntimeError("Residue numbering differs across reference tasks while deriving the embedded region")
                if occupancy_accum is None or eligible_mask is None:
                    raise RuntimeError("Embedded-region accumulator was not initialized correctly")
                occupancy_accum += np.asarray(occupancy, dtype=np.float64)
                eligible_mask = np.logical_and(eligible_mask, np.asarray(current_eligible, dtype=bool))
            source_count += 1

    if (
        reference_labels is None
        or reference_numbers is None
        or occupancy_accum is None
        or eligible_mask is None
        or half_thickness is None
        or use_curvature is None
        or curvature_radius is None
        or curvature_sign is None
    ):
        raise RuntimeError("Failed to derive the embedded residue region from the reference tasks")

    occupancy_mean = occupancy_accum / float(source_count)
    selected_mask = np.logical_and(eligible_mask, occupancy_mean >= float(occupancy_min))
    selected_indices = [int(index) for index in np.flatnonzero(selected_mask)]
    if not selected_indices:
        raise RuntimeError(
            "No residues met the embedded-region occupancy criterion; lower HYBRID_SWEEP_EMBEDDED_OCCUPANCY_MIN or inspect the reference trajectories"
        )

    return {
        "selection_mode": "reference_membrane_depth_occupancy",
        "selection_source": geometry_path,
        "half_thickness": float(half_thickness),
        "use_curvature": int(use_curvature),
        "curvature_radius": float(curvature_radius),
        "curvature_sign": float(curvature_sign),
        "occupancy_min": float(occupancy_min),
        "source_reference_task_count": int(source_count),
        "embedded_residue_indices": selected_indices,
        "embedded_residue_numbers_1based": [int(reference_numbers[index]) for index in selected_indices],
        "embedded_residue_labels": [str(reference_labels[index]) for index in selected_indices],
        "embedded_residue_occupancy": [float(occupancy_mean[index]) for index in selected_indices],
        "all_residue_numbers_1based": [int(x) for x in reference_numbers],
        "all_residue_labels": [str(x) for x in reference_labels],
        "all_residue_occupancy": [float(x) for x in occupancy_mean],
    }


def _safe_profile_correlation(a: np.ndarray, b: np.ndarray) -> float:
    if a.shape != b.shape:
        raise ValueError("Profile shapes must match for correlation")
    a_std = float(np.std(a, ddof=0))
    b_std = float(np.std(b, ddof=0))
    if a_std <= 0.0 or b_std <= 0.0:
        return 1.0 if np.allclose(a, b) else float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def _mean_std_max(values: np.ndarray) -> Dict[str, float]:
    arr = np.asarray(values, dtype=np.float64)
    return {
        "mean": float(np.mean(arr)),
        "std": float(np.std(arr, ddof=0)),
        "max": float(np.max(arr)),
    }


def _framewise_ca_radius_of_gyration(ca_xyz: np.ndarray) -> np.ndarray:
    centered = ca_xyz - np.mean(ca_xyz, axis=1, keepdims=True)
    return np.sqrt(np.mean(np.sum(centered * centered, axis=2), axis=1))


def _framewise_ca_span(ca_xyz: np.ndarray) -> np.ndarray:
    delta = ca_xyz[:, :, None, :] - ca_xyz[:, None, :, :]
    dist = np.sqrt(np.sum(delta * delta, axis=3))
    return np.max(dist, axis=(1, 2))


def _safe_ratio(value: Any, reference: Any) -> float | None:
    try:
        value = float(value)
        reference = float(reference)
    except (TypeError, ValueError):
        return None
    if not (math.isfinite(value) and math.isfinite(reference) and reference > 0.0):
        return None
    return float(value / reference)


def _format_ratio(value: float | None) -> str:
    if value is None:
        return ""
    return f"{float(value):.10g}"


def _analysis_settings_with_defaults(settings: Dict[str, Any] | None) -> Dict[str, Any]:
    merged = dict(settings or {})
    merged.setdefault("burn_in_fraction", DEFAULT_BURN_IN_FRACTION)
    merged.setdefault("trendline_samples", DEFAULT_TRENDLINE_SAMPLES)
    merged.setdefault("embedded_occupancy_min", DEFAULT_EMBEDDED_OCCUPANCY_MIN)
    merged.setdefault("rmsf_atom_role", "CA")
    merged.setdefault("alignment_iterations", 3)
    merged.setdefault(
        "stability_mean_rmsf_ratio_max",
        DEFAULT_STABILITY_MEAN_RMSF_RATIO_MAX,
    )
    merged.setdefault(
        "stability_max_rmsf_ratio_max",
        DEFAULT_STABILITY_MAX_RMSF_RATIO_MAX,
    )
    merged.setdefault(
        "stability_ca_rg_ratio_max",
        DEFAULT_STABILITY_CA_RG_RATIO_MAX,
    )
    merged.setdefault(
        "stability_ca_span_ratio_max",
        DEFAULT_STABILITY_CA_SPAN_RATIO_MAX,
    )
    return merged


def _normalize_analysis_for_assembly(analysis: Dict[str, Any]) -> Dict[str, Any]:
    normalized = dict(analysis)
    profile = np.asarray(normalized["residue_rmsf_angstrom"], dtype=np.float64)
    normalized.setdefault("max_rmsf_angstrom", float(np.max(profile)))
    normalized.setdefault("residue_numbers_1based", list(range(1, int(profile.shape[0]) + 1)))
    normalized.setdefault("n_frames_dropped_nonfinite", 0)
    normalized.setdefault("selected_backbone_atom_count", 0)
    normalized.setdefault("selected_particle_classes", [])
    normalized.setdefault("selected_atom_roles", [])
    normalized.setdefault("selected_residue_indices", [])
    normalized.setdefault("selected_residue_numbers_1based", normalized["residue_numbers_1based"])
    normalized.setdefault("selected_residue_labels", normalized.get("residue_labels", []))
    normalized.setdefault("selected_region_kind", "protein_backbone")
    normalized.setdefault("selected_region_source", "")
    normalized.setdefault("mean_embedded_region_rmsd_angstrom", None)
    normalized.setdefault("std_embedded_region_rmsd_angstrom", None)
    normalized.setdefault("max_embedded_region_rmsd_angstrom", None)
    for key in (
        "ca_radius_of_gyration_angstrom_mean",
        "ca_radius_of_gyration_angstrom_std",
        "ca_radius_of_gyration_angstrom_max",
        "ca_span_angstrom_mean",
        "ca_span_angstrom_std",
        "ca_span_angstrom_max",
    ):
        normalized.setdefault(key, None)
    return normalized


def _mean_available(values: Iterable[Any]) -> float | None:
    usable: List[float] = []
    for value in values:
        try:
            parsed = float(value)
        except (TypeError, ValueError):
            continue
        if math.isfinite(parsed):
            usable.append(parsed)
    if not usable:
        return None
    return float(np.mean(np.asarray(usable, dtype=np.float64)))


def _csv_optional_float(value: Any) -> str | float:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return ""
    if not math.isfinite(parsed):
        return ""
    return float(parsed)


def _stability_failure_summary(
    analysis: Dict[str, Any],
    reference_baselines: Dict[str, float | None],
    analysis_settings: Dict[str, Any],
) -> Dict[str, Any]:
    checks = [
        (
            "mean_rmsf_ratio",
            analysis.get("mean_rmsf_angstrom"),
            reference_baselines.get("mean_rmsf_angstrom"),
            float(analysis_settings["stability_mean_rmsf_ratio_max"]),
            "mean_rmsf",
        ),
        (
            "max_rmsf_ratio",
            analysis.get("max_rmsf_angstrom"),
            reference_baselines.get("max_rmsf_angstrom"),
            float(analysis_settings["stability_max_rmsf_ratio_max"]),
            "max_rmsf",
        ),
        (
            "ca_radius_of_gyration_ratio",
            analysis.get("ca_radius_of_gyration_angstrom_mean"),
            reference_baselines.get("ca_radius_of_gyration_angstrom_mean"),
            float(analysis_settings["stability_ca_rg_ratio_max"]),
            "ca_rg",
        ),
        (
            "ca_span_ratio",
            analysis.get("ca_span_angstrom_mean"),
            reference_baselines.get("ca_span_angstrom_mean"),
            float(analysis_settings["stability_ca_span_ratio_max"]),
            "ca_span",
        ),
    ]
    ratios: Dict[str, float | None] = {}
    reasons: List[str] = []
    for ratio_key, value, baseline, limit, label in checks:
        ratio = _safe_ratio(value, baseline)
        ratios[ratio_key] = ratio
        if ratio is not None and ratio > limit:
            reasons.append(f"{label}={ratio:.3f}>{limit:.3f}")
    return {
        "is_stable": not reasons,
        "reasons": reasons,
        "ratios": ratios,
    }


def _analyze_rmsf_profile(
    stage_file: Path,
    *,
    burn_in_fraction: float,
    alignment_iterations: int,
    embedded_residue_indices: Sequence[int] | None,
) -> Dict[str, Any]:
    with h5py.File(stage_file, "r") as h5:
        inp = h5["input"]
        if "output/pos" not in h5:
            raise ValueError(f"Missing output/pos in {stage_file}")
        pos = _normalize_output_positions(h5["output/pos"][:])
        n_frame_total = int(pos.shape[0])
        if n_frame_total < 5:
            raise ValueError(f"Need at least 5 output frames for RMSF analysis, found {n_frame_total}")

        backbone = _extract_backbone_map(h5)
        atom_index_full = np.asarray(backbone["atom_index"], dtype=np.int64)
        residue_labels_full = [str(x) for x in backbone["residue_labels"]]
        residue_numbers_full = [int(index + 1) for index in range(len(residue_labels_full))]
        selected_residue_indices = list(range(atom_index_full.shape[0]))
        selected_region_kind = "protein_backbone"
        selected_region_source = str(backbone["source"])
        if embedded_residue_indices is not None:
            selected_residue_indices = sorted({int(index) for index in embedded_residue_indices})
            if not selected_residue_indices:
                raise ValueError("Embedded residue selection is empty")
            if min(selected_residue_indices) < 0 or max(selected_residue_indices) >= atom_index_full.shape[0]:
                raise ValueError("Embedded residue selection is out of range for the residue backbone map")
            atom_index = np.asarray(atom_index_full[selected_residue_indices, :], dtype=np.int64)
            selected_region_kind = "reference_embedded_region"
            selected_region_source = "analysis_manifest.embedded_region"
        else:
            atom_index = np.asarray(atom_index_full, dtype=np.int64)
        residue_labels = [residue_labels_full[index] for index in selected_residue_indices]
        residue_numbers = [residue_numbers_full[index] for index in selected_residue_indices]

        flat_index = atom_index.reshape(-1)
        if np.max(flat_index) >= pos.shape[1]:
            raise ValueError("Backbone map refers to atoms outside output/pos")

        selected_particle_classes: List[str] = []
        if "particle_class" in inp:
            selected_particle_classes = _unique_selected_strings(inp["particle_class"][:], flat_index)
            if any(not value.startswith("PROTEIN") for value in selected_particle_classes):
                raise ValueError(
                    "RMSF backbone selection is not protein-only: "
                    + ",".join(selected_particle_classes)
                )

        selected_atom_roles: List[str] = []
        if "atom_roles" in inp:
            selected_atom_roles = _unique_selected_strings(inp["atom_roles"][:], flat_index)

        burn_start = int(math.floor(float(burn_in_fraction) * n_frame_total))
        burn_start = min(max(0, burn_start), n_frame_total - 1)
        pos = pos[burn_start:]
        n_frame_used = int(pos.shape[0])
        if n_frame_used < 5:
            raise ValueError("Too few frames remain after burn-in for RMSF analysis")

        backbone_xyz = np.asarray(pos[:, flat_index, :], dtype=np.float64)
        box = _load_box_frames_optional(h5, n_frame_total)
        box_unwrap_applied = False
        if box is not None:
            box = box[burn_start:]
            backbone_xyz = _unwrap_xyz(backbone_xyz, box)
            box_unwrap_applied = True

        frame_finite_mask = np.all(np.isfinite(backbone_xyz), axis=(1, 2))
        dropped_nonfinite_frames = int(np.count_nonzero(~frame_finite_mask))
        if dropped_nonfinite_frames:
            backbone_xyz = backbone_xyz[frame_finite_mask]
            n_frame_used = int(backbone_xyz.shape[0])
        if n_frame_used < 5:
            raise ValueError(
                f"Too few finite backbone frames remain after filtering for RMSF analysis: {n_frame_used}"
            )

        aligned_backbone, mean_backbone = _iterative_align_frames(
            backbone_xyz,
            n_iter=alignment_iterations,
        )
        n_residue = int(atom_index.shape[0])
        n_backbone_atom = int(atom_index.shape[1])
        aligned_backbone = aligned_backbone.reshape(n_frame_used, n_residue, n_backbone_atom, 3)
        ca_xyz = aligned_backbone[:, :, int(backbone["ca_column"]), :]
        mean_ca = np.mean(ca_xyz, axis=0)
        rmsf = np.sqrt(np.mean(np.sum((ca_xyz - mean_ca) ** 2, axis=2), axis=0))
        embedded_region_rmsd = np.sqrt(
            np.mean(np.sum((ca_xyz - mean_ca[None, :, :]) ** 2, axis=2), axis=1)
        )
        embedded_region_rmsd_stats = _mean_std_max(embedded_region_rmsd)
        ca_radius_of_gyration = _framewise_ca_radius_of_gyration(ca_xyz)
        ca_span = _framewise_ca_span(ca_xyz)
        ca_rg_stats = _mean_std_max(ca_radius_of_gyration)
        ca_span_stats = _mean_std_max(ca_span)

        return {
            "n_frames_total": n_frame_total,
            "n_frames_used": n_frame_used,
            "burn_in_frames": int(burn_start),
            "n_frames_dropped_nonfinite": dropped_nonfinite_frames,
            "n_residue": n_residue,
            "selected_backbone_atom_count": int(flat_index.size),
            "selected_particle_classes": selected_particle_classes,
            "selected_atom_roles": selected_atom_roles,
            "selected_region_kind": selected_region_kind,
            "selected_region_source": selected_region_source,
            "selected_residue_indices": [int(x) for x in selected_residue_indices],
            "selected_residue_numbers_1based": [int(x) for x in residue_numbers],
            "selected_residue_labels": [str(x) for x in residue_labels],
            "residue_labels": [str(x) for x in residue_labels],
            "residue_numbers_1based": [int(x) for x in residue_numbers],
            "backbone_map_source": str(backbone["source"]),
            "alignment_iterations": int(alignment_iterations),
            "box_unwrap_applied": int(box_unwrap_applied),
            "mean_rmsf_angstrom": float(np.mean(rmsf)),
            "std_rmsf_angstrom": float(np.std(rmsf, ddof=0)),
            "max_rmsf_angstrom": float(np.max(rmsf)),
            "mean_embedded_region_rmsd_angstrom": float(embedded_region_rmsd_stats["mean"]),
            "std_embedded_region_rmsd_angstrom": float(embedded_region_rmsd_stats["std"]),
            "max_embedded_region_rmsd_angstrom": float(embedded_region_rmsd_stats["max"]),
            "ca_radius_of_gyration_angstrom_mean": float(ca_rg_stats["mean"]),
            "ca_radius_of_gyration_angstrom_std": float(ca_rg_stats["std"]),
            "ca_radius_of_gyration_angstrom_max": float(ca_rg_stats["max"]),
            "ca_span_angstrom_mean": float(ca_span_stats["mean"]),
            "ca_span_angstrom_std": float(ca_span_stats["std"]),
            "ca_span_angstrom_max": float(ca_span_stats["max"]),
            "residue_rmsf_angstrom": [float(x) for x in rmsf],
            "backbone_mean_structure_shape": list(mean_backbone.shape),
        }


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
                "analysis_target_file": str(result["analysis_target_file"]),
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
    analysis_settings = _analysis_settings_with_defaults(
        manifest["settings"].get("analysis_settings")
    )
    embedded_region = _derive_reference_embedded_region(
        tasks,
        burn_in_fraction=float(analysis_settings["burn_in_fraction"]),
        occupancy_min=float(analysis_settings["embedded_occupancy_min"]),
    )
    payload = {
        "schema": SCHEMA_ANALYSIS_MANIFEST,
        "created_at_utc": _now_utc(),
        "base_dir": str(base_dir),
        "manifest_path": str(_manifest_path(base_dir)),
        "pdb_id": manifest["settings"]["pdb_id"],
        "analysis_settings": analysis_settings,
        "embedded_region": embedded_region,
        "n_analysis_targets": len(tasks),
        "tasks": tasks,
    }
    _write_json(_analysis_manifest_path(base_dir), payload)
    print(f"Analysis manifest: {_analysis_manifest_path(base_dir)}")
    print(f"Discovered analysis targets: {len(tasks)}")
    return 0


def _run_analysis_task(base_dir: Path, analysis_manifest: Dict[str, Any], task_entry: Dict[str, Any], overwrite: bool) -> Dict[str, Any]:
    task = dict(task_entry["task"])
    code = str(task["code"])
    result_path = _analysis_result_path(base_dir, code)
    if result_path.exists() and not overwrite:
        existing = _load_json(result_path)
        if existing.get("schema") == SCHEMA_ANALYSIS_RESULT and existing.get("success"):
            print(f"Skipping completed analysis task {code}")
            return existing

    stage_file = Path(task_entry["analysis_target_file"]).expanduser().resolve()
    if not stage_file.exists():
        raise RuntimeError(f"Analysis target file not found: {stage_file}")

    analysis_settings = _analysis_settings_with_defaults(
        analysis_manifest.get("analysis_settings")
    )
    embedded_region = dict(analysis_manifest.get("embedded_region") or {})
    try:
        analysis = _analyze_rmsf_profile(
            stage_file,
            burn_in_fraction=float(analysis_settings["burn_in_fraction"]),
            alignment_iterations=int(analysis_settings["alignment_iterations"]),
            embedded_residue_indices=embedded_region.get("embedded_residue_indices"),
        )
        result = {
            "schema": SCHEMA_ANALYSIS_RESULT,
            "created_at_utc": _now_utc(),
            "base_dir": str(base_dir),
            "analysis_task_id": int(task_entry["analysis_task_id"]),
            "success": True,
            "task": task,
            "analysis_target_file": str(stage_file),
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
            "analysis_target_file": str(stage_file),
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
            _run_analysis_task(base_dir, analysis_manifest, task_entry, overwrite=bool(args.overwrite))
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
    _run_analysis_task(base_dir, analysis_manifest, tasks[task_id], overwrite=bool(args.overwrite))
    return 0


def _fit_trendline(x: np.ndarray, y: np.ndarray, samples: int) -> Dict[str, Any]:
    if x.size != y.size or x.size < 1:
        raise ValueError("Trend-line fit requires at least one sampled point")
    degree = min(2, int(x.size) - 1)
    coeffs = np.polyfit(x, y, degree) if degree > 0 else np.asarray([float(np.mean(y))], dtype=np.float64)
    fitted = np.polyval(coeffs, x)
    y_mean = float(np.mean(y))
    ss_res = float(np.sum((y - fitted) ** 2))
    ss_tot = float(np.sum((y - y_mean) ** 2))
    r2 = 1.0 if ss_tot <= 0.0 else 1.0 - ss_res / ss_tot
    grid_x = np.linspace(float(np.min(x)), float(np.max(x)), int(samples))
    grid_y = np.polyval(coeffs, grid_x)
    best_index = int(np.argmin(grid_y))
    return {
        "degree": int(degree),
        "coefficients": [float(v) for v in np.asarray(coeffs, dtype=np.float64)],
        "r2": float(r2),
        "grid_x": [float(v) for v in grid_x],
        "grid_y": [float(v) for v in grid_y],
        "recommended_interface_scale": float(grid_x[best_index]),
        "recommended_metric_value": float(grid_y[best_index]),
    }


def _prepare_matplotlib(assembled_dir: Path):
    mpl_config_dir = assembled_dir / ".mplconfig"
    cache_dir = assembled_dir / ".cache"
    mpl_config_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config_dir))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def _write_interface_scale_rmsf_plot(
    *,
    assembled_dir: Path,
    condition_rows: List[Dict[str, Any]],
    stable_condition_rows: List[Dict[str, Any]],
    trend: Dict[str, Any] | None,
    recommendation: Dict[str, Any],
) -> Dict[str, Any]:
    plt = _prepare_matplotlib(assembled_dir)

    png_path = assembled_dir / "interface_scale_vs_rmsf_difference.png"
    svg_path = assembled_dir / "interface_scale_vs_rmsf_difference.svg"

    fig, ax = plt.subplots(figsize=(7.5, 4.8))

    stable_by_scale = {
        float(row["interface_scale"]): row
        for row in condition_rows
        if int(row.get("n_replicates_stable", 0)) > 0
        and str(row.get("condition_embedded_region_rmsd_delta_vs_reference_angstrom", "")) != ""
    }

    if stable_condition_rows:
        x_values = np.asarray(
            [float(row["interface_scale"]) for row in stable_condition_rows],
            dtype=np.float64,
        )
        y_values = np.asarray(
            [
                float(row["condition_embedded_region_rmsd_delta_vs_reference_angstrom"])
                for row in stable_condition_rows
            ],
            dtype=np.float64,
        )
        y_errors = np.asarray(
            [
                float(row["task_embedded_region_rmsd_delta_std_angstrom"])
                if str(row.get("task_embedded_region_rmsd_delta_std_angstrom", "")) != ""
                else 0.0
                for row in stable_condition_rows
            ],
            dtype=np.float64,
        )
        ax.errorbar(
            x_values,
            y_values,
            yerr=y_errors,
            fmt="o",
            color="#1f6f8b",
            ecolor="#7db7c7",
            elinewidth=1.2,
            capsize=3.0,
            markersize=6.5,
            label="Stable condition mean",
            zorder=3,
        )

        for row in stable_condition_rows:
            scale = float(row["interface_scale"])
            stable_meta = stable_by_scale[scale]
            label = f"{int(stable_meta['n_replicates_stable'])}/{int(stable_meta['n_replicates_completed'])}"
            ax.annotate(
                label,
                (scale, float(row["condition_embedded_region_rmsd_delta_vs_reference_angstrom"])),
                textcoords="offset points",
                xytext=(0, 6),
                ha="center",
                fontsize=8,
                color="#24505c",
            )

    excluded_scales = [
        float(row["interface_scale"])
        for row in condition_rows
        if int(row.get("condition_excluded_from_fit", 0)) == 1
    ]
    if excluded_scales:
        excluded_y = 1.0
        if stable_condition_rows:
            excluded_y = float(
                max(
                    float(row["condition_embedded_region_rmsd_delta_vs_reference_angstrom"])
                    for row in stable_condition_rows
                )
                * 1.06
            )
        ax.scatter(
            excluded_scales,
            np.full(len(excluded_scales), excluded_y, dtype=np.float64),
            marker="x",
            s=48,
            linewidths=1.6,
            color="#c84b31",
            label="Excluded from fit",
            zorder=4,
        )
        for scale in excluded_scales:
            excluded_meta = next(row for row in condition_rows if float(row["interface_scale"]) == scale)
            label = f"0/{int(excluded_meta['n_replicates_completed'])}"
            ax.annotate(
                label,
                (scale, excluded_y),
                textcoords="offset points",
                xytext=(0, 6),
                ha="center",
                fontsize=8,
                color="#8c2d19",
            )

    if trend is not None:
        grid_x = np.asarray(trend["grid_x"], dtype=np.float64)
        grid_y = np.asarray(trend["grid_y"], dtype=np.float64)
        ax.plot(
            grid_x,
            grid_y,
            color="#d95f02",
            linewidth=2.0,
            label=f"Quadratic fit (R^2={float(trend['r2']):.3f})" if int(trend["degree"]) == 2 else "Trend fit",
            zorder=2,
        )
        ax.axvline(
            float(trend["recommended_interface_scale"]),
            color="#d95f02",
            linestyle="--",
            linewidth=1.2,
            alpha=0.8,
        )

    best_sampled_scale = recommendation.get("best_sampled_interface_scale")
    best_sampled_rmse = recommendation.get(
        "best_sampled_condition_embedded_region_rmsd_delta_vs_reference_angstrom"
    )
    if best_sampled_scale is not None and best_sampled_rmse is not None:
        ax.scatter(
            [float(best_sampled_scale)],
            [float(best_sampled_rmse)],
            marker="*",
            s=140,
            color="#2a9d8f",
            edgecolors="#1d5c54",
            linewidths=0.8,
            label="Best sampled",
            zorder=5,
        )

    ax.set_xlabel("PROTEIN_ENV_INTERFACE_SCALE")
    ax.set_ylabel("Embedded-Region RMSD Delta Vs Reference (Angstrom)")
    ax.set_title("1rkl Embedded-Region RMSD Calibration")
    ax.grid(True, alpha=0.25, linewidth=0.7)
    if condition_rows:
        scales = sorted(float(row["interface_scale"]) for row in condition_rows)
        ax.set_xlim(min(scales) - 0.03, max(scales) + 0.03)
    ax.legend(frameon=False, fontsize=8, loc="best")

    note = "Labels show stable/completed replicate counts per scale"
    ax.text(
        0.99,
        0.02,
        note,
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=8,
        color="#4b5563",
    )

    fig.tight_layout()
    fig.savefig(png_path, dpi=200)
    fig.savefig(svg_path)
    plt.close(fig)
    return {
        "png": str(png_path),
        "svg": str(svg_path),
    }


def _select_profile_comparison_scale(
    *,
    stable_scales: Iterable[float],
    best_sampled_scale: Any,
    trendline_scale: Any,
) -> Dict[str, Any]:
    sampled = sorted({float(value) for value in stable_scales})
    if not sampled:
        return {
            "selected_interface_scale": None,
            "target_interface_scale": None,
            "selection_basis": "",
        }

    target_scale: float | None = None
    selection_basis = ""
    if trendline_scale is not None:
        target_scale = float(trendline_scale)
        selection_basis = "trendline_recommendation"
    elif best_sampled_scale is not None:
        target_scale = float(best_sampled_scale)
        selection_basis = "best_sampled"

    if target_scale is None:
        target_scale = float(sampled[0])
        selection_basis = "first_stable_sample"

    selected_scale = min(
        sampled,
        key=lambda value: (abs(float(value) - target_scale), abs(float(value)), float(value)),
    )
    if selection_basis == "trendline_recommendation" and not math.isclose(
        float(selected_scale), float(target_scale), rel_tol=0.0, abs_tol=1e-12
    ):
        selection_basis = "nearest_stable_sample_to_trendline_recommendation"

    return {
        "selected_interface_scale": float(selected_scale),
        "target_interface_scale": float(target_scale),
        "selection_basis": selection_basis,
    }


def _write_best_scale_reference_rmsf_plot(
    *,
    assembled_dir: Path,
    reference_numbers: Sequence[int],
    reference_labels: Sequence[str],
    reference_mean: np.ndarray,
    reference_std: np.ndarray,
    profile_meta: Dict[str, Any],
    selection: Dict[str, Any],
) -> Dict[str, Any]:
    if selection.get("selected_interface_scale") is None:
        return {}

    selected_scale = float(selection["selected_interface_scale"])
    n_stable = int(profile_meta["n_replicates_stable"])
    n_completed = int(profile_meta["n_replicates_completed"])
    selection_note = str(selection.get("selection_basis", "")).replace("_", " ")
    note = (
        f"Profile scale selection: {selection_note}; "
        f"stable/completed replicates = {n_stable}/{n_completed}"
        if selection_note
        else f"Stable/completed replicates = {n_stable}/{n_completed}"
    )
    return _write_reference_comparison_plot(
        assembled_dir=assembled_dir,
        png_path=assembled_dir / "best_interface_scale_rmsf_vs_reference.png",
        svg_path=assembled_dir / "best_interface_scale_rmsf_vs_reference.svg",
        reference_numbers=reference_numbers,
        reference_labels=reference_labels,
        reference_mean=reference_mean,
        reference_std=reference_std,
        condition_mean=np.asarray(profile_meta["mean"], dtype=np.float64),
        condition_std=np.asarray(profile_meta["std"], dtype=np.float64),
        interface_scale=selected_scale,
        note=note,
    )


def _write_reference_comparison_plot(
    *,
    assembled_dir: Path,
    png_path: Path,
    svg_path: Path,
    reference_numbers: Sequence[int],
    reference_labels: Sequence[str],
    reference_mean: np.ndarray,
    reference_std: np.ndarray,
    condition_mean: np.ndarray,
    condition_std: np.ndarray,
    interface_scale: float,
    note: str,
) -> Dict[str, Any]:
    plt = _prepare_matplotlib(assembled_dir)

    x_values = np.asarray([int(x) for x in reference_numbers], dtype=np.int64)
    fig, ax = plt.subplots(figsize=(8.4, 4.8))

    ax.fill_between(
        x_values,
        reference_mean - reference_std,
        reference_mean + reference_std,
        color="#b8c4d6",
        alpha=0.45,
        linewidth=0.0,
        label="Reference ±1σ",
    )
    ax.plot(
        x_values,
        reference_mean,
        color="#3d5a80",
        linewidth=2.1,
        label="Reference mean",
    )

    ax.fill_between(
        x_values,
        condition_mean - condition_std,
        condition_mean + condition_std,
        color="#f4a261",
        alpha=0.35,
        linewidth=0.0,
        label=f"Scale {interface_scale:.3f} ±1σ",
    )
    ax.plot(
        x_values,
        condition_mean,
        color="#d1495b",
        linewidth=2.1,
        label=f"Scale {interface_scale:.3f} mean",
    )

    tick_positions = [int(x_values[0])]
    for position in x_values[1:]:
        if (int(position) - int(tick_positions[-1])) >= 5:
            tick_positions.append(int(position))
    if int(x_values[-1]) not in tick_positions:
        tick_positions.append(int(x_values[-1]))
    tick_positions = sorted(set(tick_positions))
    tick_label_map = {
        int(number): f"{int(number)}\n{str(label)}"
        for number, label in zip(reference_numbers, reference_labels)
    }
    tick_labels = [tick_label_map.get(int(pos), str(int(pos))) for pos in tick_positions]

    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, fontsize=8)
    ax.set_xlabel("Residue")
    ax.set_ylabel("Backbone RMSF (Angstrom)")
    ax.set_title(f"1rkl Embedded-Region RMSF: Reference vs Interface Scale {interface_scale:.3f}")
    ax.grid(True, alpha=0.22, linewidth=0.7)
    ax.legend(frameon=False, fontsize=8, loc="upper right")

    if note:
        ax.text(
            0.99,
            0.02,
            note,
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=8,
            color="#4b5563",
        )

    fig.tight_layout()
    fig.savefig(png_path, dpi=200)
    fig.savefig(svg_path)
    plt.close(fig)
    return {
        "png": str(png_path),
        "svg": str(svg_path),
    }


def _write_all_scale_reference_rmsf_plots(
    *,
    assembled_dir: Path,
    reference_numbers: Sequence[int],
    reference_labels: Sequence[str],
    reference_mean: np.ndarray,
    reference_std: np.ndarray,
    condition_rows: Sequence[Dict[str, Any]],
    stable_profile_by_scale: Dict[float, Dict[str, Any]],
) -> Dict[str, Any]:
    plot_dir = assembled_dir / "scale_rmsf_vs_reference"
    plot_dir.mkdir(parents=True, exist_ok=True)
    index_csv = assembled_dir / "scale_rmsf_vs_reference_index.csv"

    rows: List[Dict[str, Any]] = []
    for row in sorted(condition_rows, key=lambda item: float(item["interface_scale"])):
        interface_scale = float(row["interface_scale"])
        profile_meta = stable_profile_by_scale.get(interface_scale)
        plot_png = ""
        plot_svg = ""
        plot_error = ""
        if profile_meta is not None:
            tag = _format_float_tag(interface_scale, digits=4)
            png_path = plot_dir / f"scale{tag}_rmsf_vs_reference.png"
            svg_path = plot_dir / f"scale{tag}_rmsf_vs_reference.svg"
            note = (
                f"Sampled scale overlay; stable/completed replicates = "
                f"{int(profile_meta['n_replicates_stable'])}/{int(profile_meta['n_replicates_completed'])}"
            )
            try:
                outputs = _write_reference_comparison_plot(
                    assembled_dir=assembled_dir,
                    png_path=png_path,
                    svg_path=svg_path,
                    reference_numbers=reference_numbers,
                    reference_labels=reference_labels,
                    reference_mean=reference_mean,
                    reference_std=reference_std,
                    condition_mean=np.asarray(profile_meta["mean"], dtype=np.float64),
                    condition_std=np.asarray(profile_meta["std"], dtype=np.float64),
                    interface_scale=interface_scale,
                    note=note,
                )
                plot_png = outputs.get("png", "")
                plot_svg = outputs.get("svg", "")
            except Exception as exc:
                plot_error = str(exc)

        rows.append(
            {
                "interface_scale": float(interface_scale),
                "n_replicates_completed": int(row["n_replicates_completed"]),
                "n_replicates_stable": int(row["n_replicates_stable"]),
                "condition_embedded_region_rmsd_delta_vs_reference_angstrom": row.get(
                    "condition_embedded_region_rmsd_delta_vs_reference_angstrom", ""
                ),
                "plot_png": plot_png,
                "plot_svg": plot_svg,
                "plot_error": plot_error,
            }
        )

    with index_csv.open("w", encoding="utf-8", newline="") as fh:
        fieldnames = [
            "interface_scale",
            "n_replicates_completed",
            "n_replicates_stable",
            "condition_embedded_region_rmsd_delta_vs_reference_angstrom",
            "plot_png",
            "plot_svg",
            "plot_error",
        ]
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    return {
        "plot_dir": str(plot_dir),
        "index_csv": str(index_csv),
    }


def assemble_analysis(base_dir: Path) -> int:
    manifest = _load_manifest(base_dir)
    analysis_manifest = _load_analysis_manifest(base_dir)
    results_dir = _analysis_results_dir(base_dir)
    assembled_dir = _analysis_assembled_dir(base_dir)
    assembled_dir.mkdir(parents=True, exist_ok=True)

    successful_results = _successful_analysis_results_from_dir(results_dir)
    for result in successful_results:
        result["analysis"] = _normalize_analysis_for_assembly(result["analysis"])
    reference_results = [item for item in successful_results if item["task"]["kind"] == "reference"]
    hybrid_results = [item for item in successful_results if item["task"]["kind"] == "hybrid"]
    if not reference_results:
        raise RuntimeError("No successful reference RMSF analyses found")
    if not hybrid_results:
        raise RuntimeError("No successful hybrid RMSF analyses found")

    reference_profiles = np.asarray(
        [result["analysis"]["residue_rmsf_angstrom"] for result in reference_results],
        dtype=np.float64,
    )
    reference_labels = [str(x) for x in reference_results[0]["analysis"]["residue_labels"]]
    reference_numbers = [
        int(x)
        for x in reference_results[0]["analysis"].get(
            "residue_numbers_1based",
            list(range(1, int(reference_profiles.shape[1]) + 1)),
        )
    ]
    n_residue = int(reference_profiles.shape[1])
    for result in successful_results:
        profile = np.asarray(result["analysis"]["residue_rmsf_angstrom"], dtype=np.float64)
        labels = [str(x) for x in result["analysis"]["residue_labels"]]
        numbers = [
            int(x)
            for x in result["analysis"].get(
                "residue_numbers_1based",
                list(range(1, int(profile.shape[0]) + 1)),
            )
        ]
        if profile.shape[0] != n_residue:
            raise RuntimeError("Residue-count mismatch across embedded-region RMSF analysis results")
        if labels != reference_labels:
            raise RuntimeError("Residue-label mismatch across embedded-region RMSF analysis results")
        if numbers != reference_numbers:
            raise RuntimeError("Residue-number mismatch across embedded-region RMSF analysis results")

    reference_mean = np.mean(reference_profiles, axis=0)
    reference_std = np.std(reference_profiles, axis=0, ddof=0)
    reference_baselines = {
        "mean_rmsf_angstrom": _mean_available(
            item["analysis"].get("mean_rmsf_angstrom") for item in reference_results
        ),
        "max_rmsf_angstrom": _mean_available(
            item["analysis"].get("max_rmsf_angstrom") for item in reference_results
        ),
        "mean_embedded_region_rmsd_angstrom": _mean_available(
            item["analysis"].get("mean_embedded_region_rmsd_angstrom") for item in reference_results
        ),
        "ca_radius_of_gyration_angstrom_mean": _mean_available(
            item["analysis"].get("ca_radius_of_gyration_angstrom_mean") for item in reference_results
        ),
        "ca_span_angstrom_mean": _mean_available(
            item["analysis"].get("ca_span_angstrom_mean") for item in reference_results
        ),
    }
    analysis_settings = _analysis_settings_with_defaults(
        analysis_manifest.get("analysis_settings")
    )
    embedded_region = dict(analysis_manifest.get("embedded_region") or {})

    reference_profile_csv = assembled_dir / "reference_profile.csv"
    with reference_profile_csv.open("w", encoding="utf-8", newline="") as fh:
        fieldnames = [
            "residue_index",
            "residue_number_1based",
            "residue_label",
            "reference_rmsf_mean_angstrom",
            "reference_rmsf_std_angstrom",
        ]
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for index, label in enumerate(reference_labels, start=1):
            writer.writerow(
                {
                    "residue_index": index,
                    "residue_number_1based": int(reference_numbers[index - 1]),
                    "residue_label": label,
                    "reference_rmsf_mean_angstrom": float(reference_mean[index - 1]),
                    "reference_rmsf_std_angstrom": float(reference_std[index - 1]),
                }
            )

    task_rows: List[Dict[str, Any]] = []
    profile_rows: List[Dict[str, Any]] = []
    grouped_all: Dict[float, List[Dict[str, Any]]] = {}
    grouped_stable: Dict[float, List[Dict[str, Any]]] = {}
    for result in successful_results:
        task = result["task"]
        analysis = result["analysis"]
        profile = np.asarray(analysis["residue_rmsf_angstrom"], dtype=np.float64)
        stability_is_stable = True
        stability_reasons: List[str] = []
        stability_ratios = {
            "mean_rmsf_ratio": None,
            "max_rmsf_ratio": None,
            "ca_radius_of_gyration_ratio": None,
            "ca_span_ratio": None,
        }
        if task["kind"] == "hybrid":
            stability = _stability_failure_summary(
                analysis=analysis,
                reference_baselines=reference_baselines,
                analysis_settings=analysis_settings,
            )
            stability_is_stable = bool(stability["is_stable"])
            stability_reasons = list(stability["reasons"])
            stability_ratios = dict(stability["ratios"])
        row: Dict[str, Any] = {
            "analysis_task_id": int(result["analysis_task_id"]),
            "task_id": int(task.get("task_id", -1)),
            "code": str(task["code"]),
            "kind": str(task["kind"]),
            "method": str(task["method"]),
            "interface_scale": (
                ""
                if task["kind"] != "hybrid"
                else f"{float(task['interface_scale']):.10g}"
            ),
            "replicate": int(task["replicate"]),
            "selected_region_kind": str(analysis.get("selected_region_kind", "")),
            "selected_region_source": str(analysis.get("selected_region_source", "")),
            "selected_residue_count": int(len(analysis.get("selected_residue_numbers_1based", []))),
            "mean_embedded_region_rmsd_angstrom": _csv_optional_float(
                analysis.get("mean_embedded_region_rmsd_angstrom")
            ),
            "std_embedded_region_rmsd_angstrom": _csv_optional_float(
                analysis.get("std_embedded_region_rmsd_angstrom")
            ),
            "max_embedded_region_rmsd_angstrom": _csv_optional_float(
                analysis.get("max_embedded_region_rmsd_angstrom")
            ),
            "mean_rmsf_angstrom": float(analysis["mean_rmsf_angstrom"]),
            "std_rmsf_angstrom": float(analysis["std_rmsf_angstrom"]),
            "max_rmsf_angstrom": float(analysis["max_rmsf_angstrom"]),
            "ca_radius_of_gyration_angstrom_mean": _csv_optional_float(
                analysis.get("ca_radius_of_gyration_angstrom_mean")
            ),
            "ca_radius_of_gyration_angstrom_std": _csv_optional_float(
                analysis.get("ca_radius_of_gyration_angstrom_std")
            ),
            "ca_radius_of_gyration_angstrom_max": _csv_optional_float(
                analysis.get("ca_radius_of_gyration_angstrom_max")
            ),
            "ca_span_angstrom_mean": _csv_optional_float(analysis.get("ca_span_angstrom_mean")),
            "ca_span_angstrom_std": _csv_optional_float(analysis.get("ca_span_angstrom_std")),
            "ca_span_angstrom_max": _csv_optional_float(analysis.get("ca_span_angstrom_max")),
            "n_frames_used": int(analysis["n_frames_used"]),
            "n_frames_dropped_nonfinite": int(analysis.get("n_frames_dropped_nonfinite", 0)),
            "selected_backbone_atom_count": int(analysis.get("selected_backbone_atom_count", 0)),
            "selected_particle_classes": ",".join(analysis.get("selected_particle_classes", [])),
            "selected_atom_roles": ",".join(analysis.get("selected_atom_roles", [])),
            "is_stable_protein_trajectory": int(stability_is_stable),
            "stability_failure_reasons": ";".join(stability_reasons),
            "mean_rmsf_ratio_to_reference": _format_ratio(stability_ratios["mean_rmsf_ratio"]),
            "max_rmsf_ratio_to_reference": _format_ratio(stability_ratios["max_rmsf_ratio"]),
            "ca_radius_of_gyration_ratio_to_reference": _format_ratio(
                stability_ratios["ca_radius_of_gyration_ratio"]
            ),
            "ca_span_ratio_to_reference": _format_ratio(stability_ratios["ca_span_ratio"]),
            "analysis_target_file": str(result["analysis_target_file"]),
            "embedded_region_rmsd_signed_delta_vs_reference_angstrom": "",
            "embedded_region_rmsd_delta_vs_reference_angstrom": "",
        }
        if task["kind"] == "hybrid":
            task_embedded_rmsd = float(analysis["mean_embedded_region_rmsd_angstrom"])
            reference_embedded_rmsd = float(reference_baselines["mean_embedded_region_rmsd_angstrom"])
            signed_delta = task_embedded_rmsd - reference_embedded_rmsd
            row["embedded_region_rmsd_signed_delta_vs_reference_angstrom"] = float(signed_delta)
            row["embedded_region_rmsd_delta_vs_reference_angstrom"] = float(abs(signed_delta))
            grouped_all.setdefault(float(task["interface_scale"]), []).append(
                {
                    "row": row,
                    "profile": profile,
                    "embedded_region_rmsd": task_embedded_rmsd,
                }
            )
            if stability_is_stable:
                grouped_stable.setdefault(float(task["interface_scale"]), []).append(
                    {
                        "row": row,
                        "profile": profile,
                        "embedded_region_rmsd": task_embedded_rmsd,
                    }
                )
        task_rows.append(row)

        for index, label in enumerate(reference_labels, start=1):
            profile_rows.append(
                {
                    "kind": str(task["kind"]),
                    "method": str(task["method"]),
                    "code": str(task["code"]),
                    "interface_scale": (
                        ""
                        if task["kind"] != "hybrid"
                        else f"{float(task['interface_scale']):.10g}"
                    ),
                    "replicate": int(task["replicate"]),
                    "is_stable_protein_trajectory": int(stability_is_stable),
                    "residue_index": index,
                    "residue_number_1based": int(reference_numbers[index - 1]),
                    "residue_label": label,
                    "rmsf_angstrom": float(profile[index - 1]),
                    "reference_rmsf_mean_angstrom": float(reference_mean[index - 1]),
                    "reference_rmsf_std_angstrom": float(reference_std[index - 1]),
                    "delta_vs_reference_angstrom": float(profile[index - 1] - reference_mean[index - 1]),
                }
            )

    task_csv = assembled_dir / "task_results.csv"
    with task_csv.open("w", encoding="utf-8", newline="") as fh:
        fieldnames = [
            "analysis_task_id",
            "task_id",
            "code",
            "kind",
            "method",
            "interface_scale",
            "replicate",
            "selected_region_kind",
            "selected_region_source",
            "selected_residue_count",
            "mean_embedded_region_rmsd_angstrom",
            "std_embedded_region_rmsd_angstrom",
            "max_embedded_region_rmsd_angstrom",
            "mean_rmsf_angstrom",
            "std_rmsf_angstrom",
            "max_rmsf_angstrom",
            "ca_radius_of_gyration_angstrom_mean",
            "ca_radius_of_gyration_angstrom_std",
            "ca_radius_of_gyration_angstrom_max",
            "ca_span_angstrom_mean",
            "ca_span_angstrom_std",
            "ca_span_angstrom_max",
            "embedded_region_rmsd_signed_delta_vs_reference_angstrom",
            "embedded_region_rmsd_delta_vs_reference_angstrom",
            "n_frames_used",
            "n_frames_dropped_nonfinite",
            "selected_backbone_atom_count",
            "selected_particle_classes",
            "selected_atom_roles",
            "is_stable_protein_trajectory",
            "stability_failure_reasons",
            "mean_rmsf_ratio_to_reference",
            "max_rmsf_ratio_to_reference",
            "ca_radius_of_gyration_ratio_to_reference",
            "ca_span_ratio_to_reference",
            "analysis_target_file",
        ]
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in sorted(task_rows, key=lambda item: (item["analysis_task_id"], item["task_id"], item["code"])):
            writer.writerow(row)

    profile_csv = assembled_dir / "residue_rmsf_profiles.csv"
    with profile_csv.open("w", encoding="utf-8", newline="") as fh:
        fieldnames = [
            "kind",
            "method",
            "code",
            "interface_scale",
            "replicate",
            "is_stable_protein_trajectory",
            "residue_index",
            "residue_number_1based",
            "residue_label",
            "rmsf_angstrom",
            "reference_rmsf_mean_angstrom",
            "reference_rmsf_std_angstrom",
            "delta_vs_reference_angstrom",
        ]
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in profile_rows:
            writer.writerow(row)

    condition_rows: List[Dict[str, Any]] = []
    condition_profile_rows: List[Dict[str, Any]] = []
    stable_condition_rows: List[Dict[str, Any]] = []
    stable_profile_by_scale: Dict[float, Dict[str, Any]] = {}
    for interface_scale, rows in sorted(grouped_all.items()):
        stable_rows = grouped_stable.get(float(interface_scale), [])
        summary: Dict[str, Any] = {
            "interface_scale": float(interface_scale),
            "n_replicates_expected": int(manifest["settings"]["hybrid_replicates"]),
            "n_replicates_completed": len(rows),
            "n_replicates_stable": len(stable_rows),
            "n_replicates_filtered_unstable": len(rows) - len(stable_rows),
            "condition_excluded_from_fit": int(len(stable_rows) == 0),
            "reference_embedded_region_rmsd_mean_angstrom": _csv_optional_float(
                reference_baselines.get("mean_embedded_region_rmsd_angstrom")
            ),
            "condition_embedded_region_rmsd_mean_angstrom": "",
            "condition_embedded_region_rmsd_std_angstrom": "",
            "condition_embedded_region_rmsd_delta_vs_reference_angstrom": "",
            "task_embedded_region_rmsd_delta_mean_angstrom": "",
            "task_embedded_region_rmsd_delta_std_angstrom": "",
        }
        if stable_rows:
            profiles = np.asarray([item["profile"] for item in stable_rows], dtype=np.float64)
            condition_mean = np.mean(profiles, axis=0)
            condition_std = np.std(profiles, axis=0, ddof=0)
            task_embedded_rmsd = np.asarray(
                [float(item["embedded_region_rmsd"]) for item in stable_rows],
                dtype=np.float64,
            )
            condition_embedded_rmsd_mean = float(np.mean(task_embedded_rmsd))
            condition_embedded_rmsd_std = float(np.std(task_embedded_rmsd, ddof=0))
            condition_embedded_rmsd_delta = abs(
                condition_embedded_rmsd_mean - float(reference_baselines["mean_embedded_region_rmsd_angstrom"])
            )
            task_delta = np.asarray(
                [float(item["row"]["embedded_region_rmsd_delta_vs_reference_angstrom"]) for item in stable_rows],
                dtype=np.float64,
            )
            summary.update(
                {
                    "condition_embedded_region_rmsd_mean_angstrom": float(condition_embedded_rmsd_mean),
                    "condition_embedded_region_rmsd_std_angstrom": float(condition_embedded_rmsd_std),
                    "condition_embedded_region_rmsd_delta_vs_reference_angstrom": float(
                        condition_embedded_rmsd_delta
                    ),
                    "task_embedded_region_rmsd_delta_mean_angstrom": float(np.mean(task_delta)),
                    "task_embedded_region_rmsd_delta_std_angstrom": float(np.std(task_delta, ddof=0)),
                }
            )
            stable_condition_rows.append(summary)
            stable_profile_by_scale[float(interface_scale)] = {
                "mean": np.asarray(condition_mean, dtype=np.float64),
                "std": np.asarray(condition_std, dtype=np.float64),
                "n_replicates_stable": int(len(stable_rows)),
                "n_replicates_completed": int(len(rows)),
                "condition_embedded_region_rmsd_delta_vs_reference_angstrom": float(
                    summary["condition_embedded_region_rmsd_delta_vs_reference_angstrom"]
                ),
            }

            for index, label in enumerate(reference_labels, start=1):
                condition_profile_rows.append(
                    {
                        "interface_scale": float(interface_scale),
                        "residue_index": index,
                        "residue_number_1based": int(reference_numbers[index - 1]),
                        "residue_label": label,
                        "condition_rmsf_mean_angstrom": float(condition_mean[index - 1]),
                        "condition_rmsf_std_angstrom": float(condition_std[index - 1]),
                        "reference_rmsf_mean_angstrom": float(reference_mean[index - 1]),
                        "reference_rmsf_std_angstrom": float(reference_std[index - 1]),
                        "delta_vs_reference_angstrom": float(
                            condition_mean[index - 1] - reference_mean[index - 1]
                        ),
                    }
                )
        condition_rows.append(summary)

    summary_csv = assembled_dir / "condition_summary.csv"
    with summary_csv.open("w", encoding="utf-8", newline="") as fh:
        fieldnames = [
            "interface_scale",
            "n_replicates_expected",
            "n_replicates_completed",
            "n_replicates_stable",
            "n_replicates_filtered_unstable",
            "condition_excluded_from_fit",
            "reference_embedded_region_rmsd_mean_angstrom",
            "condition_embedded_region_rmsd_mean_angstrom",
            "condition_embedded_region_rmsd_std_angstrom",
            "condition_embedded_region_rmsd_delta_vs_reference_angstrom",
            "task_embedded_region_rmsd_delta_mean_angstrom",
            "task_embedded_region_rmsd_delta_std_angstrom",
        ]
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in condition_rows:
            writer.writerow(row)

    condition_profile_csv = assembled_dir / "condition_profiles.csv"
    with condition_profile_csv.open("w", encoding="utf-8", newline="") as fh:
        fieldnames = [
            "interface_scale",
            "residue_index",
            "residue_number_1based",
            "residue_label",
            "condition_rmsf_mean_angstrom",
            "condition_rmsf_std_angstrom",
            "reference_rmsf_mean_angstrom",
            "reference_rmsf_std_angstrom",
            "delta_vs_reference_angstrom",
        ]
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in condition_profile_rows:
            writer.writerow(row)

    trend_csv = assembled_dir / "trendline_points.csv"
    trend: Dict[str, Any] | None = None
    with trend_csv.open("w", encoding="utf-8", newline="") as fh:
        fieldnames = [
            "interface_scale",
            "fitted_condition_embedded_region_rmsd_delta_vs_reference_angstrom",
        ]
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        if stable_condition_rows:
            trend = _fit_trendline(
                np.asarray([float(row["interface_scale"]) for row in stable_condition_rows], dtype=np.float64),
                np.asarray(
                    [
                        float(row["condition_embedded_region_rmsd_delta_vs_reference_angstrom"])
                        for row in stable_condition_rows
                    ],
                    dtype=np.float64,
                ),
                int(analysis_settings["trendline_samples"]),
            )
            for x_value, y_value in zip(trend["grid_x"], trend["grid_y"]):
                writer.writerow(
                    {
                        "interface_scale": float(x_value),
                        "fitted_condition_embedded_region_rmsd_delta_vs_reference_angstrom": float(y_value),
                    }
                )

    failed_rows: List[Dict[str, Any]] = []
    for path in sorted(results_dir.glob("*.json")):
        payload = _load_json(path)
        if payload.get("success"):
            continue
        task = payload.get("task", {})
        failed_rows.append(
            {
                "code": str(task.get("code", path.stem)),
                "kind": str(task.get("kind", "")),
                "interface_scale": (
                    ""
                    if task.get("kind") != "hybrid"
                    else f"{float(task.get('interface_scale')):.10g}"
                ),
                "replicate": int(task.get("replicate", 0) or 0),
                "failure_stage": "analysis_error",
                "error": str(payload.get("error", "")),
            }
        )
    for row in task_rows:
        if row["kind"] != "hybrid" or int(row["is_stable_protein_trajectory"]) == 1:
            continue
        failed_rows.append(
            {
                "code": str(row["code"]),
                "kind": str(row["kind"]),
                "interface_scale": str(row["interface_scale"]),
                "replicate": int(row["replicate"]),
                "failure_stage": "protein_stability_filter",
                "error": str(row["stability_failure_reasons"]),
            }
        )

    failed_csv = assembled_dir / "failed_tasks.csv"
    with failed_csv.open("w", encoding="utf-8", newline="") as fh:
        fieldnames = ["code", "kind", "interface_scale", "replicate", "failure_stage", "error"]
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in failed_rows:
            writer.writerow(row)

    best_sampled: Dict[str, Any] | None = None
    if stable_condition_rows:
        best_sampled = min(
            stable_condition_rows,
            key=lambda item: (
                float(item["condition_embedded_region_rmsd_delta_vs_reference_angstrom"]),
                float(item["interface_scale"]),
            ),
        )
    analysis_status = "ok" if stable_condition_rows else "no_stable_hybrid_trajectories"
    recommendation = {
        "created_at_utc": _now_utc(),
        "analysis_status": analysis_status,
        "metric": "condition_embedded_region_rmsd_delta_vs_reference_angstrom",
        "reference_task_count": int(len(reference_results)),
        "hybrid_condition_count": int(len(condition_rows)),
        "hybrid_condition_count_used_for_fit": int(len(stable_condition_rows)),
        "hybrid_task_count_filtered_unstable": int(
            sum(1 for row in task_rows if row["kind"] == "hybrid" and int(row["is_stable_protein_trajectory"]) == 0)
        ),
        "embedded_region": embedded_region,
        "best_sampled_interface_scale": (
            None if best_sampled is None else float(best_sampled["interface_scale"])
        ),
        "best_sampled_condition_embedded_region_rmsd_delta_vs_reference_angstrom": (
            None
            if best_sampled is None
            else float(best_sampled["condition_embedded_region_rmsd_delta_vs_reference_angstrom"])
        ),
        "trendline_recommended_interface_scale": (
            None if trend is None else float(trend["recommended_interface_scale"])
        ),
        "trendline_recommended_condition_embedded_region_rmsd_delta_vs_reference_angstrom": (
            None if trend is None else float(trend["recommended_metric_value"])
        ),
        "trendline_fit_degree": None if trend is None else int(trend["degree"]),
        "trendline_fit_coefficients": [] if trend is None else [float(x) for x in trend["coefficients"]],
        "trendline_fit_r2": None if trend is None else float(trend["r2"]),
        "note": (
            ""
            if stable_condition_rows
            else "No stable hybrid trajectories remained after protein stability filtering"
        ),
        "reference_stability_baselines": reference_baselines,
        "stability_thresholds": {
            "stability_mean_rmsf_ratio_max": float(analysis_settings["stability_mean_rmsf_ratio_max"]),
            "stability_max_rmsf_ratio_max": float(analysis_settings["stability_max_rmsf_ratio_max"]),
            "stability_ca_rg_ratio_max": float(analysis_settings["stability_ca_rg_ratio_max"]),
            "stability_ca_span_ratio_max": float(analysis_settings["stability_ca_span_ratio_max"]),
        },
    }
    profile_selection = _select_profile_comparison_scale(
        stable_scales=[float(row["interface_scale"]) for row in stable_condition_rows],
        best_sampled_scale=recommendation["best_sampled_interface_scale"],
        trendline_scale=recommendation["trendline_recommended_interface_scale"],
    )
    plot_outputs: Dict[str, str] = {}
    all_scale_plot_outputs: Dict[str, str] = {}
    best_profile_plot_outputs: Dict[str, str] = {}
    plot_error = ""
    all_scale_plot_error = ""
    best_profile_plot_error = ""
    try:
        plot_outputs = _write_interface_scale_rmsf_plot(
            assembled_dir=assembled_dir,
            condition_rows=condition_rows,
            stable_condition_rows=stable_condition_rows,
            trend=trend,
            recommendation=recommendation,
        )
    except Exception as exc:
        plot_error = str(exc)
    try:
        all_scale_plot_outputs = _write_all_scale_reference_rmsf_plots(
            assembled_dir=assembled_dir,
            reference_numbers=reference_numbers,
            reference_labels=reference_labels,
            reference_mean=np.asarray(reference_mean, dtype=np.float64),
            reference_std=np.asarray(reference_std, dtype=np.float64),
            condition_rows=condition_rows,
            stable_profile_by_scale=stable_profile_by_scale,
        )
    except Exception as exc:
        all_scale_plot_error = str(exc)
    if profile_selection.get("selected_interface_scale") is not None:
        try:
            best_profile_plot_outputs = _write_best_scale_reference_rmsf_plot(
                assembled_dir=assembled_dir,
                reference_numbers=reference_numbers,
                reference_labels=reference_labels,
                reference_mean=np.asarray(reference_mean, dtype=np.float64),
                reference_std=np.asarray(reference_std, dtype=np.float64),
                profile_meta=stable_profile_by_scale[float(profile_selection["selected_interface_scale"])],
                selection=profile_selection,
            )
        except Exception as exc:
            best_profile_plot_error = str(exc)
    recommendation["profile_comparison_target_interface_scale"] = profile_selection.get(
        "target_interface_scale"
    )
    recommendation["profile_comparison_interface_scale"] = profile_selection.get(
        "selected_interface_scale"
    )
    recommendation["profile_comparison_selection_basis"] = profile_selection.get(
        "selection_basis", ""
    )
    recommendation["scale_rmsf_vs_reference_dir"] = all_scale_plot_outputs.get("plot_dir", "")
    recommendation["scale_rmsf_vs_reference_index_csv"] = all_scale_plot_outputs.get("index_csv", "")
    recommendation["interface_scale_vs_rmsf_difference_png"] = plot_outputs.get("png", "")
    recommendation["interface_scale_vs_rmsf_difference_svg"] = plot_outputs.get("svg", "")
    recommendation["best_interface_scale_rmsf_vs_reference_png"] = best_profile_plot_outputs.get("png", "")
    recommendation["best_interface_scale_rmsf_vs_reference_svg"] = best_profile_plot_outputs.get("svg", "")
    recommendation["plot_error"] = plot_error
    recommendation["all_scale_plot_error"] = all_scale_plot_error
    recommendation["best_profile_plot_error"] = best_profile_plot_error
    _write_json(assembled_dir / "recommendation_summary.json", recommendation)

    payload = {
        "created_at_utc": _now_utc(),
        "analysis_status": analysis_status,
        "base_dir": str(base_dir),
        "pdb_id": manifest["settings"]["pdb_id"],
        "metric": recommendation["metric"],
        "embedded_region": embedded_region,
        "reference_embedded_region_rmsd_mean_angstrom": reference_baselines.get(
            "mean_embedded_region_rmsd_angstrom"
        ),
        "n_analysis_tasks_total": int(len(analysis_manifest["tasks"])),
        "n_analysis_tasks_completed_successfully": int(len(successful_results)),
        "n_reference_tasks_completed_successfully": int(len(reference_results)),
        "n_hybrid_tasks_completed_successfully": int(len(hybrid_results)),
        "n_hybrid_tasks_filtered_unstable": int(
            sum(1 for row in task_rows if row["kind"] == "hybrid" and int(row["is_stable_protein_trajectory"]) == 0)
        ),
        "n_hybrid_tasks_stable": int(
            sum(1 for row in task_rows if row["kind"] == "hybrid" and int(row["is_stable_protein_trajectory"]) == 1)
        ),
        "n_hybrid_conditions_completed": int(len(condition_rows)),
        "n_hybrid_conditions_used_for_fit": int(len(stable_condition_rows)),
        "trendline_available": int(bool(stable_condition_rows)),
        "task_results_csv": str(task_csv),
        "residue_rmsf_profiles_csv": str(profile_csv),
        "reference_profile_csv": str(reference_profile_csv),
        "condition_profiles_csv": str(condition_profile_csv),
        "condition_summary_csv": str(summary_csv),
        "trendline_points_csv": str(trend_csv),
        "failed_tasks_csv": str(failed_csv),
        "scale_rmsf_vs_reference_dir": all_scale_plot_outputs.get("plot_dir", ""),
        "scale_rmsf_vs_reference_index_csv": all_scale_plot_outputs.get("index_csv", ""),
        "interface_scale_vs_rmsf_difference_png": plot_outputs.get("png", ""),
        "interface_scale_vs_rmsf_difference_svg": plot_outputs.get("svg", ""),
        "best_interface_scale_rmsf_vs_reference_png": best_profile_plot_outputs.get("png", ""),
        "best_interface_scale_rmsf_vs_reference_svg": best_profile_plot_outputs.get("svg", ""),
        "profile_comparison_interface_scale": profile_selection.get("selected_interface_scale"),
        "profile_comparison_target_interface_scale": profile_selection.get("target_interface_scale"),
        "profile_comparison_selection_basis": profile_selection.get("selection_basis", ""),
        "plot_error": plot_error,
        "all_scale_plot_error": all_scale_plot_error,
        "best_profile_plot_error": best_profile_plot_error,
    }
    _write_json(assembled_dir / "summary.json", payload)
    print(f"Assembled analysis results written under {assembled_dir}")
    return 0


def cmd_assemble_analysis(args: argparse.Namespace) -> int:
    base_dir = Path(args.base_dir).expanduser().resolve()
    return assemble_analysis(base_dir)


def _array_script_content(round_manifest_path: Path, base_dir: Path, n_tasks: int) -> str:
    train_walltime = os.environ.get("HYBRID_SWEEP_TRAIN_WALLTIME", "24:00:00").strip() or "24:00:00"
    cpus_per_task = int(os.environ.get("HYBRID_SWEEP_CPUS_PER_TASK", "1"))
    directives = "\n".join(_optional_sbatch_directives("TRAIN"))
    if directives:
        directives += "\n"
    return f"""#!/bin/bash
#SBATCH --job-name=hybrid_rmsf_sweep
#SBATCH --output={_slurm_dir(base_dir) / "array_%A_%a.out"}
#SBATCH --error={_slurm_dir(base_dir) / "array_%A_%a.err"}
#SBATCH --array=0-{max(0, n_tasks - 1)}
#SBATCH --time={train_walltime}
#SBATCH --cpus-per-task={cpus_per_task}
{directives}set -euo pipefail

PROJECT_ROOT="{REPO_ROOT}"
WORKFLOW_DIR="{WORKFLOW_DIR}"
ROUND_MANIFEST="{round_manifest_path}"
TASK_ID="${{SLURM_ARRAY_TASK_ID:?}}"

if [ -f /etc/profile.d/modules.sh ]; then
  source /etc/profile.d/modules.sh
fi

if command -v module >/dev/null 2>&1; then
  module load python/3.11.9 || true
  module load cmake || true
  module load openmpi || true
  module load "${{HYBRID_SWEEP_HDF5_MODULE:-hdf5/1.14.3}}" || true
fi

if [ -f "$PROJECT_ROOT/.venv/bin/activate" ]; then
  source "$PROJECT_ROOT/.venv/bin/activate"
fi

PYTHON_BIN="${{HYBRID_SWEEP_PYTHON:-python3}}"
export HYBRID_SWEEP_PROJECT_ROOT="$PROJECT_ROOT"
export PYTHONUNBUFFERED=1

"$PYTHON_BIN" -u "$WORKFLOW_DIR/workflow.py" run-array-task --round-manifest "$ROUND_MANIFEST" --task-id "$TASK_ID"
"""


def _collector_script_content(base_dir: Path) -> str:
    collect_walltime = os.environ.get("HYBRID_SWEEP_COLLECT_WALLTIME", "01:00:00").strip() or "01:00:00"
    directives = "\n".join(_optional_sbatch_directives("COLLECT"))
    if directives:
        directives += "\n"
    return f"""#!/bin/bash
#SBATCH --job-name=hybrid_rmsf_collect
#SBATCH --output={_slurm_dir(base_dir) / "collect_%j.out"}
#SBATCH --error={_slurm_dir(base_dir) / "collect_%j.err"}
#SBATCH --time={collect_walltime}
#SBATCH --cpus-per-task=1
{directives}set -euo pipefail

PROJECT_ROOT="{REPO_ROOT}"
WORKFLOW_DIR="{WORKFLOW_DIR}"
BASE_DIR="{base_dir}"

if [ -f /etc/profile.d/modules.sh ]; then
  source /etc/profile.d/modules.sh
fi

if command -v module >/dev/null 2>&1; then
  module load python/3.11.9 || true
  module load cmake || true
  module load openmpi || true
  module load "${{HYBRID_SWEEP_HDF5_MODULE:-hdf5/1.14.3}}" || true
fi

if [ -f "$PROJECT_ROOT/.venv/bin/activate" ]; then
  source "$PROJECT_ROOT/.venv/bin/activate"
fi

PYTHON_BIN="${{HYBRID_SWEEP_PYTHON:-python3}}"
export HYBRID_SWEEP_PROJECT_ROOT="$PROJECT_ROOT"
export PYTHONUNBUFFERED=1

"$PYTHON_BIN" -u "$WORKFLOW_DIR/workflow.py" assemble-results --base-dir "$BASE_DIR"
"""


def _analysis_array_script_content(round_manifest_path: Path, base_dir: Path, n_tasks: int) -> str:
    walltime = os.environ.get("HYBRID_SWEEP_ANALYSIS_WALLTIME", "04:00:00").strip() or "04:00:00"
    cpus_per_task = int(os.environ.get("HYBRID_SWEEP_ANALYSIS_CPUS_PER_TASK", "1"))
    directives = "\n".join(_optional_sbatch_directives("ANALYSIS"))
    if directives:
        directives += "\n"
    return f"""#!/bin/bash
#SBATCH --job-name=hybrid_rmsf_analysis
#SBATCH --output={_analysis_slurm_dir(base_dir) / "array_%A_%a.out"}
#SBATCH --error={_analysis_slurm_dir(base_dir) / "array_%A_%a.err"}
#SBATCH --array=0-{max(0, n_tasks - 1)}
#SBATCH --time={walltime}
#SBATCH --cpus-per-task={cpus_per_task}
{directives}set -euo pipefail

PROJECT_ROOT="{REPO_ROOT}"
WORKFLOW_DIR="{WORKFLOW_DIR}"
ROUND_MANIFEST="{round_manifest_path}"
TASK_ID="${{SLURM_ARRAY_TASK_ID:?}}"

if [ -f /etc/profile.d/modules.sh ]; then
  source /etc/profile.d/modules.sh
fi

if command -v module >/dev/null 2>&1; then
  module load python/3.11.9 || true
  module load cmake || true
  module load openmpi || true
  module load "${{HYBRID_SWEEP_HDF5_MODULE:-hdf5/1.14.3}}" || true
fi

if [ -f "$PROJECT_ROOT/.venv/bin/activate" ]; then
  source "$PROJECT_ROOT/.venv/bin/activate"
fi

PYTHON_BIN="${{HYBRID_SWEEP_PYTHON:-python3}}"
export HYBRID_SWEEP_PROJECT_ROOT="$PROJECT_ROOT"
export PYTHONUNBUFFERED=1

"$PYTHON_BIN" -u "$WORKFLOW_DIR/workflow.py" run-analysis-task --round-manifest "$ROUND_MANIFEST" --task-id "$TASK_ID"
"""


def _analysis_collector_script_content(base_dir: Path) -> str:
    collect_walltime = os.environ.get("HYBRID_SWEEP_ANALYSIS_COLLECT_WALLTIME", "01:00:00").strip() or "01:00:00"
    directives = "\n".join(_optional_sbatch_directives("ANALYSIS_COLLECT"))
    if directives:
        directives += "\n"
    return f"""#!/bin/bash
#SBATCH --job-name=hybrid_rmsf_analysis_collect
#SBATCH --output={_analysis_slurm_dir(base_dir) / "collect_%j.out"}
#SBATCH --error={_analysis_slurm_dir(base_dir) / "collect_%j.err"}
#SBATCH --time={collect_walltime}
#SBATCH --cpus-per-task=1
{directives}set -euo pipefail

PROJECT_ROOT="{REPO_ROOT}"
WORKFLOW_DIR="{WORKFLOW_DIR}"
BASE_DIR="{base_dir}"

if [ -f /etc/profile.d/modules.sh ]; then
  source /etc/profile.d/modules.sh
fi

if command -v module >/dev/null 2>&1; then
  module load python/3.11.9 || true
  module load cmake || true
  module load openmpi || true
  module load "${{HYBRID_SWEEP_HDF5_MODULE:-hdf5/1.14.3}}" || true
fi

if [ -f "$PROJECT_ROOT/.venv/bin/activate" ]; then
  source "$PROJECT_ROOT/.venv/bin/activate"
fi

PYTHON_BIN="${{HYBRID_SWEEP_PYTHON:-python3}}"
export HYBRID_SWEEP_PROJECT_ROOT="$PROJECT_ROOT"
export PYTHONUNBUFFERED=1

"$PYTHON_BIN" -u "$WORKFLOW_DIR/workflow.py" assemble-analysis --base-dir "$BASE_DIR"
"""


def cmd_submit_slurm(args: argparse.Namespace) -> int:
    base_dir = Path(args.base_dir).expanduser().resolve()
    manifest = _load_manifest(base_dir)
    tasks = manifest["tasks"]
    if not tasks:
        raise RuntimeError("No tasks available to submit")

    slurm_dir = _slurm_dir(base_dir)
    slurm_dir.mkdir(parents=True, exist_ok=True)
    round_manifest_path = slurm_dir / "round_manifest.json"
    array_script = slurm_dir / "run_array.sbatch"
    collect_script = slurm_dir / "collect_results.sbatch"

    round_manifest = {
        "schema": SCHEMA_SLURM_ROUND,
        "created_at_utc": _now_utc(),
        "base_dir": str(base_dir),
        "manifest_path": str(_manifest_path(base_dir)),
        "task_count": int(len(tasks)),
    }
    _write_json(round_manifest_path, round_manifest)
    _write_text(array_script, _array_script_content(round_manifest_path, base_dir, len(tasks)), executable=True)
    _write_text(collect_script, _collector_script_content(base_dir), executable=True)

    if args.no_submit:
        print(f"Staged Slurm run scripts under {slurm_dir}")
        return 0

    sbatch = _shutil_which("sbatch")
    if sbatch is None:
        raise RuntimeError("sbatch not found in PATH")
    array_job_id = _run_sbatch([sbatch, "--parsable", str(array_script)])
    collect_job_id = _run_sbatch([sbatch, "--parsable", f"--dependency=afterok:{array_job_id}", str(collect_script)])
    _write_json(
        slurm_dir / "submission.json",
        {
            "created_at_utc": _now_utc(),
            "array_job_id": array_job_id,
            "collect_job_id": collect_job_id,
        },
    )
    print(f"Submitted array job {array_job_id}")
    print(f"Submitted collector job {collect_job_id}")
    return 0


def cmd_submit_analysis_slurm(args: argparse.Namespace) -> int:
    base_dir = Path(args.base_dir).expanduser().resolve()
    analysis_manifest = _load_analysis_manifest(base_dir)
    tasks = analysis_manifest["tasks"]
    if not tasks:
        raise RuntimeError("No analysis tasks available to submit")

    slurm_dir = _analysis_slurm_dir(base_dir)
    slurm_dir.mkdir(parents=True, exist_ok=True)
    round_manifest_path = slurm_dir / "round_manifest.json"
    array_script = slurm_dir / "analyze_array.sbatch"
    collect_script = slurm_dir / "collect_analysis.sbatch"

    round_manifest = {
        "schema": SCHEMA_ANALYSIS_SLURM,
        "created_at_utc": _now_utc(),
        "base_dir": str(base_dir),
        "analysis_manifest_path": str(_analysis_manifest_path(base_dir)),
        "task_count": int(len(tasks)),
    }
    _write_json(round_manifest_path, round_manifest)
    _write_text(
        array_script,
        _analysis_array_script_content(round_manifest_path, base_dir, len(tasks)),
        executable=True,
    )
    _write_text(collect_script, _analysis_collector_script_content(base_dir), executable=True)

    if args.no_submit:
        print(f"Staged analysis Slurm scripts under {slurm_dir}")
        return 0

    sbatch = _shutil_which("sbatch")
    if sbatch is None:
        raise RuntimeError("sbatch not found in PATH")
    array_job_id = _run_sbatch([sbatch, "--parsable", str(array_script)])
    collect_job_id = _run_sbatch([sbatch, "--parsable", f"--dependency=afterok:{array_job_id}", str(collect_script)])
    _write_json(
        slurm_dir / "submission.json",
        {
            "created_at_utc": _now_utc(),
            "array_job_id": array_job_id,
            "collect_job_id": collect_job_id,
        },
    )
    print(f"Submitted analysis array job {array_job_id}")
    print(f"Submitted analysis collector job {collect_job_id}")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Hybrid interface RMSF calibration sweep")
    sub = parser.add_subparsers(dest="command", required=True)

    p_init = sub.add_parser("init-run", help="Initialize a new interface-scale RMSF sweep")
    p_init.add_argument("--base-dir", default=str(WORKFLOW_DIR / "runs" / "default"))
    p_init.add_argument("--pdb-id", default=DEFAULT_PDB_ID)
    p_init.add_argument("--interface-scales", type=_float_list_arg, default=DEFAULT_INTERFACE_SCALES)
    p_init.add_argument("--hybrid-replicates", type=int, default=DEFAULT_HYBRID_REPLICATES)
    p_init.add_argument("--reference-replicates", type=int, default=DEFAULT_REFERENCE_REPLICATES)
    p_init.add_argument("--seed", type=int, default=DEFAULT_SEED)
    p_init.add_argument("--burn-in-fraction", type=float, default=DEFAULT_BURN_IN_FRACTION)
    p_init.add_argument("--trendline-samples", type=int, default=DEFAULT_TRENDLINE_SAMPLES)
    p_init.set_defaults(func=cmd_init_run)

    p_local = sub.add_parser("run-local", help="Run sweep tasks locally")
    p_local.add_argument("--base-dir", default=str(WORKFLOW_DIR / "runs" / "default"))
    p_local.add_argument("--max-tasks", type=int)
    p_local.add_argument("--start-task", type=int, default=0)
    p_local.add_argument("--overwrite", action="store_true")
    p_local.add_argument("--no-assemble", action="store_true")
    p_local.set_defaults(func=cmd_run_local)

    p_array = sub.add_parser("run-array-task", help="Run one task from a Slurm round manifest")
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

    p_init_analysis = sub.add_parser("init-analysis", help="Discover completed outputs for RMSF analysis")
    p_init_analysis.add_argument("--base-dir", default=str(WORKFLOW_DIR / "runs" / "default"))
    p_init_analysis.set_defaults(func=cmd_init_analysis)

    p_analysis_local = sub.add_parser("run-analysis-local", help="Run completed analysis tasks locally")
    p_analysis_local.add_argument("--base-dir", default=str(WORKFLOW_DIR / "runs" / "default"))
    p_analysis_local.add_argument("--max-tasks", type=int)
    p_analysis_local.add_argument("--start-task", type=int, default=0)
    p_analysis_local.add_argument("--overwrite", action="store_true")
    p_analysis_local.add_argument("--no-assemble", action="store_true")
    p_analysis_local.set_defaults(func=cmd_run_analysis_local)

    p_analysis_task = sub.add_parser("run-analysis-task", help="Run one analysis task from a Slurm round manifest")
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
