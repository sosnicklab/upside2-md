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
import sys
import tarfile
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Sequence, Tuple

import h5py
import numpy as np


WORKFLOW_DIR = Path(os.path.abspath(__file__)).parent
REPO_ROOT = WORKFLOW_DIR.parent
PREP_SCRIPT = REPO_ROOT / "py" / "martini_prepare_system.py"
UPSIDE_EXECUTABLE = REPO_ROOT / "obj" / "upside"
DEFAULT_BILAYER_PDB = REPO_ROOT / "parameters" / "dryMARTINI" / "DOPC.pdb"

SCHEMA_MANIFEST = "bilayer_lateral_diffusion_manifest_v1"
SCHEMA_SLURM_ROUND = "bilayer_lateral_diffusion_slurm_round_v1"
SCHEMA_TASK_RESULT = "bilayer_lateral_diffusion_task_result_v1"
SCHEMA_STAGE7_ANALYSIS_MANIFEST = "bilayer_lateral_diffusion_stage7_analysis_manifest_v1"
SCHEMA_STAGE7_ANALYSIS_SLURM = "bilayer_lateral_diffusion_stage7_analysis_slurm_v1"
SCHEMA_STAGE7_ANALYSIS_RESULT = "bilayer_lateral_diffusion_stage7_analysis_result_v1"

DEFAULT_DAMPING_VALUES = [4.0, 6.0, 8.0, 12.0, 16.0, 32.0, 64.0]
DEFAULT_MASS_SCALES = [1.0, 0.75, 0.5, 0.25, 0.1]
DEFAULT_TEMPERATURE_VALUES = [0.70, 0.80, 0.90, 1.00, 1.10]

DEFAULT_REPLICATES = 3
DEFAULT_EQ_STEPS = 10000
DEFAULT_PROD_STEPS = 100000
DEFAULT_MIN_60_MAX_ITER = 500
DEFAULT_MIN_61_MAX_ITER = 500
DEFAULT_TIME_STEP = 0.01
DEFAULT_INTEGRATOR = "v"
DEFAULT_MAX_FORCE = 0.0
DEFAULT_EQ_FRAME_STEPS = 1000
DEFAULT_PROD_FRAME_STEPS = 100
DEFAULT_SEED = 20260406
DEFAULT_XY_SCALE = 1.0
DEFAULT_SALT_MOLAR = 0.15
DEFAULT_ION_CUTOFF = 4.0
DEFAULT_BOX_PADDING_Z = 20.0

UPSIDE_MARTINI_ENERGY_CONVERSION = 2.914952774272
UPSIDE_MARTINI_LENGTH_CONVERSION = 10.0
UPSIDE_COMPRESSIBILITY = 14.521180763676
UPSIDE_BAR_1_TO_EUP_PER_A3 = 0.000020659477
UPSIDE_NPT_TAU = 4.0
UPSIDE_NPT_INTERVAL = 10
UPSIDE_NPT_SEMI = 1
UPSIDE_NPT_DEBUG = 1
THERMOSTAT_INTERVAL = -1
BASE_PDB_ID = "bilayer_diffusion"


def _now_utc() -> str:
    return datetime.now(timezone.utc).isoformat()


def _manifest_path(base_dir: Path) -> Path:
    return base_dir / "diffusion_manifest.json"


def _slurm_dir(base_dir: Path) -> Path:
    return base_dir / "slurm"


def _system_dir(base_dir: Path) -> Path:
    return base_dir / "system"


def _task_dir(base_dir: Path, code: str) -> Path:
    return base_dir / "tasks" / code


def _results_dir(base_dir: Path) -> Path:
    return base_dir / "results"


def _task_results_dir(base_dir: Path) -> Path:
    return _results_dir(base_dir) / "tasks"


def _task_result_path(base_dir: Path, code: str) -> Path:
    return _task_results_dir(base_dir) / f"{code}.json"


def _assembled_dir(base_dir: Path) -> Path:
    return base_dir / "assembled"


def _stage7_analysis_dir(base_dir: Path) -> Path:
    return base_dir / "stage7-analysis"


def _stage7_analysis_manifest_path(base_dir: Path) -> Path:
    return _stage7_analysis_dir(base_dir) / "analysis_manifest.json"


def _stage7_analysis_results_dir(base_dir: Path) -> Path:
    return _stage7_analysis_dir(base_dir) / "results" / "tasks"


def _stage7_analysis_result_path(base_dir: Path, code: str) -> Path:
    return _stage7_analysis_results_dir(base_dir) / f"{code}.json"


def _stage7_analysis_assembled_dir(base_dir: Path) -> Path:
    return _stage7_analysis_dir(base_dir) / "assembled"


def _stage7_analysis_slurm_dir(base_dir: Path) -> Path:
    return _stage7_analysis_dir(base_dir) / "slurm"


def _download_dir(base_dir: Path) -> Path:
    return base_dir / "download"


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


def _load_stage7_analysis_manifest(base_dir: Path) -> Dict[str, Any]:
    path = _stage7_analysis_manifest_path(base_dir)
    if not path.exists():
        raise RuntimeError(f"Stage-7 analysis manifest not found: {path}. Run init-stage7-analysis first.")
    manifest = _load_json(path)
    if manifest.get("schema") != SCHEMA_STAGE7_ANALYSIS_MANIFEST:
        raise RuntimeError(f"Unexpected stage-7 analysis manifest schema in {path}: {manifest.get('schema')}")
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
            f"BILAYER_DIFF_{phase}_{key}",
            os.environ.get(f"BILAYER_DIFF_SBATCH_{key}", ""),
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
    env["UPSIDE_MARTINI_ENERGY_CONVERSION"] = str(UPSIDE_MARTINI_ENERGY_CONVERSION)
    env["UPSIDE_MARTINI_LENGTH_CONVERSION"] = str(UPSIDE_MARTINI_LENGTH_CONVERSION)
    env["UPSIDE_NPT_TARGET_PXY"] = str(UPSIDE_BAR_1_TO_EUP_PER_A3)
    env["UPSIDE_NPT_TARGET_PZ"] = str(UPSIDE_BAR_1_TO_EUP_PER_A3)
    env["UPSIDE_NPT_TAU"] = str(UPSIDE_NPT_TAU)
    env["UPSIDE_NPT_COMPRESSIBILITY"] = str(UPSIDE_COMPRESSIBILITY)
    env["UPSIDE_NPT_COMPRESSIBILITY_XY"] = str(UPSIDE_COMPRESSIBILITY)
    env["UPSIDE_NPT_COMPRESSIBILITY_Z"] = "0.0"
    env["UPSIDE_NPT_INTERVAL"] = str(UPSIDE_NPT_INTERVAL)
    env["UPSIDE_NPT_SEMI"] = str(UPSIDE_NPT_SEMI)
    env["UPSIDE_NPT_DEBUG"] = str(UPSIDE_NPT_DEBUG)
    env["UPSIDE_EWALD_ENABLE"] = "1"
    env["UPSIDE_EWALD_ALPHA"] = "0.2"
    env["UPSIDE_EWALD_KMAX"] = "5"
    env["UPSIDE_MARTINI_FF_DIR"] = str(REPO_ROOT / "parameters" / "dryMARTINI")
    env["PYTHONUNBUFFERED"] = "1"
    _prepend_env_path(env, "PATH", str(REPO_ROOT / "obj"))
    _prepend_env_path(env, "PYTHONPATH", str(REPO_ROOT / "py"))
    return env


def _require_runtime_files() -> None:
    if not PREP_SCRIPT.exists():
        raise RuntimeError(f"Preparation script not found: {PREP_SCRIPT}")
    if not UPSIDE_EXECUTABLE.exists():
        raise RuntimeError(f"Upside executable not found: {UPSIDE_EXECUTABLE}")
    if not DEFAULT_BILAYER_PDB.exists():
        raise RuntimeError(f"Default DOPC template not found: {DEFAULT_BILAYER_PDB}")


def _stage_schedule(manifest: Dict[str, Any]) -> List[Dict[str, Any]]:
    settings = manifest["settings"]
    eq_steps = int(settings["equilibration_steps"])
    prod_steps = int(settings["production_steps"])
    return [
        {
            "label": "6.0",
            "prepare_stage": "minimization",
            "run_kind": "minimization",
            "max_iter": int(settings["min_60_max_iter"]),
            "lipidhead_fc": 0.0,
            "steps": 0,
        },
        {
            "label": "6.1",
            "prepare_stage": "npt_prod",
            "run_kind": "minimization",
            "max_iter": int(settings["min_61_max_iter"]),
            "lipidhead_fc": 0.0,
            "steps": 0,
        },
        {
            "label": "6.2",
            "prepare_stage": "npt_equil",
            "run_kind": "md",
            "lipidhead_fc": 200.0,
            "steps": eq_steps,
            "frame_steps": int(settings["equilibration_frame_steps"]),
        },
        {
            "label": "6.3",
            "prepare_stage": "npt_equil_reduced",
            "run_kind": "md",
            "lipidhead_fc": 100.0,
            "steps": eq_steps,
            "frame_steps": int(settings["equilibration_frame_steps"]),
        },
        {
            "label": "6.4",
            "prepare_stage": "npt_prod",
            "run_kind": "md",
            "lipidhead_fc": 50.0,
            "steps": eq_steps,
            "frame_steps": int(settings["equilibration_frame_steps"]),
        },
        {
            "label": "6.5",
            "prepare_stage": "npt_prod",
            "run_kind": "md",
            "lipidhead_fc": 20.0,
            "steps": eq_steps,
            "frame_steps": int(settings["equilibration_frame_steps"]),
        },
        {
            "label": "6.6",
            "prepare_stage": "npt_prod",
            "run_kind": "md",
            "lipidhead_fc": 10.0,
            "steps": eq_steps,
            "frame_steps": int(settings["equilibration_frame_steps"]),
        },
        {
            "label": "7.0",
            "prepare_stage": "npt_prod",
            "run_kind": "md",
            "lipidhead_fc": 0.0,
            "steps": prod_steps,
            "frame_steps": int(settings["production_frame_steps"]),
        },
    ]


def _python_bin() -> str:
    return sys.executable or "python3"


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


def _prepare_shared_bilayer(manifest: Dict[str, Any]) -> Dict[str, Any]:
    base_dir = Path(manifest["base_dir"]).expanduser().resolve()
    system_dir = _system_dir(base_dir)
    runtime_pdb = system_dir / "dopc_bilayer.MARTINI.pdb"
    summary_json = system_dir / "prep_summary.json"
    prep_run_dir = system_dir / "prep_run"
    system_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        _python_bin(),
        str(PREP_SCRIPT),
        "--mode",
        "bilayer",
        "--pdb-id",
        BASE_PDB_ID,
        "--runtime-pdb-output",
        str(runtime_pdb),
        "--prepare-structure",
        "1",
        "--run-dir",
        str(prep_run_dir),
        "--summary-json",
        str(summary_json),
        "--bilayer-pdb",
        str(Path(manifest["settings"]["bilayer_pdb"]).expanduser().resolve()),
        "--xy-scale",
        str(manifest["settings"]["xy_scale"]),
        "--box-padding-z",
        str(manifest["settings"]["box_padding_z"]),
        "--salt-molar",
        str(manifest["settings"]["salt_molar"]),
        "--ion-cutoff",
        str(manifest["settings"]["ion_cutoff"]),
        "--seed",
        str(manifest["settings"]["seed"]),
    ]
    env = _base_runtime_env()
    log_path = system_dir / "prepare.log"
    _call_logged(cmd, log_path, env)
    summary = _load_json(summary_json) if summary_json.exists() else {}
    return {
        "runtime_pdb": str(runtime_pdb),
        "prep_summary_json": str(summary_json),
        "prepare_log": str(log_path),
        "summary": summary,
    }


def _build_manifest(args: argparse.Namespace, base_dir: Path) -> Dict[str, Any]:
    damping_values = sorted(dict.fromkeys(float(x) for x in args.damping_values))
    mass_scales = sorted(dict.fromkeys(float(x) for x in args.mass_scales), reverse=True)
    temperature_values = sorted(dict.fromkeys(float(x) for x in args.temperature_values))
    replicates = int(args.replicates)
    settings = {
        "bilayer_pdb": str(Path(args.bilayer_pdb).expanduser().resolve()),
        "seed": int(args.seed),
        "xy_scale": float(args.xy_scale),
        "box_padding_z": float(args.box_padding_z),
        "salt_molar": float(args.salt_molar),
        "ion_cutoff": float(args.ion_cutoff),
        "damping_values": damping_values,
        "mass_scales": mass_scales,
        "temperature_values": temperature_values,
        "replicates": replicates,
        "min_60_max_iter": int(args.min_60_max_iter),
        "min_61_max_iter": int(args.min_61_max_iter),
        "equilibration_steps": int(args.equilibration_steps),
        "production_steps": int(args.production_steps),
        "time_step": float(args.time_step),
        "integrator": str(args.integrator),
        "max_force": float(args.max_force),
        "equilibration_frame_steps": int(args.equilibration_frame_steps),
        "production_frame_steps": int(args.production_frame_steps),
    }
    tasks: List[Dict[str, Any]] = []
    task_id = 0
    for temperature in temperature_values:
        for damping in damping_values:
            for mass_scale in mass_scales:
                for replicate in range(replicates):
                    code = (
                        f"tau{_format_float_tag(damping)}_"
                        f"m{_format_float_tag(mass_scale)}_"
                        f"t{_format_float_tag(temperature)}_"
                        f"r{replicate + 1:02d}"
                    )
                    tasks.append(
                        {
                            "task_id": task_id,
                            "code": code,
                            "temperature": float(temperature),
                            "damping": float(damping),
                            "mass_scale": float(mass_scale),
                            "replicate": int(replicate + 1),
                            "seed": _task_seed(int(args.seed), task_id),
                        }
                    )
                    task_id += 1
    return {
        "schema": SCHEMA_MANIFEST,
        "created_at_utc": _now_utc(),
        "base_dir": str(base_dir),
        "workflow_dir": str(WORKFLOW_DIR),
        "repo_root": str(REPO_ROOT),
        "pdb_id": BASE_PDB_ID,
        "settings": settings,
        "tasks": tasks,
    }


def cmd_init_run(args: argparse.Namespace) -> int:
    _require_runtime_files()
    base_dir = Path(args.base_dir).expanduser().resolve()
    base_dir.mkdir(parents=True, exist_ok=True)
    manifest = _build_manifest(args, base_dir)
    system_info = _prepare_shared_bilayer(manifest)
    manifest["system"] = system_info
    _write_json(_manifest_path(base_dir), manifest)
    print(f"Initialized bilayer diffusion run: {base_dir}")
    print(f"Shared runtime PDB: {system_info['runtime_pdb']}")
    print(f"Task count: {len(manifest['tasks'])}")
    return 0


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


def _compute_msd(xy: np.ndarray, times: np.ndarray, max_lag: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
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


def _decode_str_array(dataset: h5py.Dataset) -> np.ndarray:
    arr = dataset[:]
    out = []
    for item in arr:
        if isinstance(item, (bytes, np.bytes_)):
            out.append(item.decode("utf-8").strip())
        else:
            out.append(str(item).strip())
    return np.asarray(out, dtype=object)


def _load_box_frames(h5: h5py.File, n_frame: int) -> np.ndarray:
    if "output/box" in h5:
        return _normalize_box(h5["output/box"][:])
    grp = h5["input/potential/martini_potential"]
    box = np.array([grp.attrs["x_len"], grp.attrs["y_len"], grp.attrs["z_len"]], dtype=np.float64)
    return np.repeat(box[None, :], n_frame, axis=0)


def _wrapped_fraction(values: np.ndarray, box_length: float) -> np.ndarray:
    return np.mod(values, box_length) / box_length


def _select_lipid_atoms_and_po4(
    atom_names: np.ndarray,
    molecule_ids: np.ndarray,
    particle_class: np.ndarray | None,
) -> Tuple[np.ndarray, np.ndarray, str]:
    po4_mask = atom_names == "PO4"
    if particle_class is not None:
        lipid_mask = particle_class == "LIPID"
        lipid_po4_mask = lipid_mask & po4_mask
        if np.any(lipid_po4_mask):
            return lipid_mask, lipid_po4_mask, "particle_class"

    if not np.any(po4_mask):
        raise ValueError("No PO4 beads found in production stage file")

    lipid_molecule_ids = np.unique(molecule_ids[po4_mask])
    lipid_mask = np.isin(molecule_ids, lipid_molecule_ids)
    lipid_po4_mask = lipid_mask & po4_mask
    if not np.any(lipid_po4_mask):
        raise ValueError("Could not infer lipid molecules from PO4 beads")
    return lipid_mask, lipid_po4_mask, "po4_molecule_fallback"


def _split_leaflets_from_wrapped_z(
    initial_po4_z: np.ndarray,
    initial_box_z: float,
) -> Tuple[np.ndarray, np.ndarray, float, float]:
    z_frac = _wrapped_fraction(initial_po4_z, initial_box_z)
    order = np.argsort(z_frac)
    sorted_frac = z_frac[order]
    cyclic_gaps = np.diff(np.concatenate([sorted_frac, sorted_frac[:1] + 1.0]))
    cut_index = int(np.argmax(cyclic_gaps))
    cut_origin = float(sorted_frac[(cut_index + 1) % sorted_frac.shape[0]])

    rotated = np.mod(z_frac - cut_origin, 1.0)
    rotated_order = np.argsort(rotated)
    rotated_sorted = rotated[rotated_order]
    if rotated_sorted.shape[0] < 2:
        raise ValueError("Need at least two PO4 beads to split bilayer leaflets")

    split_gaps = np.diff(rotated_sorted)
    split_index = int(np.argmax(split_gaps))
    split_threshold = float(0.5 * (rotated_sorted[split_index] + rotated_sorted[split_index + 1]))
    bottom_mask = rotated <= split_threshold
    top_mask = ~bottom_mask
    if not np.any(top_mask) or not np.any(bottom_mask):
        raise ValueError("Could not split PO4 beads into two wrapped-z leaflets")
    return top_mask, bottom_mask, cut_origin, split_threshold


def _load_reference_positions_and_box(stage_file: Path) -> Tuple[np.ndarray, float]:
    with h5py.File(stage_file, "r") as h5:
        if "output/pos" in h5:
            pos = _normalize_output_positions(h5["output/pos"][:])[0]
        else:
            pos = np.asarray(h5["input/pos"][:], dtype=np.float64)
        if "output/box" in h5:
            box_z = float(_normalize_box(h5["output/box"][:])[0, 2])
        else:
            grp = h5["input/potential/martini_potential"]
            box_z = float(grp.attrs["z_len"])
    return pos, box_z


def _analyze_production(stage_file: Path, reference_stage_file: Path | None = None) -> Dict[str, Any]:
    with h5py.File(stage_file, "r") as h5:
        if "output/pos" not in h5:
            raise ValueError(f"Missing output/pos in {stage_file}")
        pos = _normalize_output_positions(h5["output/pos"][:])
        n_frame, n_atom, _ = pos.shape
        if n_frame < 20:
            raise ValueError(f"Need at least 20 output frames for MSD analysis, found {n_frame}")

        if "output/time" in h5:
            times = _normalize_time(h5["output/time"][:])
        else:
            times = np.arange(n_frame, dtype=np.float64)
        if times.shape[0] != n_frame:
            raise ValueError("output/time length does not match output/pos frames")

        box = _load_box_frames(h5, n_frame)
        atom_names = _decode_str_array(h5["input/atom_names"])
        molecule_ids = np.asarray(h5["input/molecule_ids"][:], dtype=np.int32)
        particle_class = None
        if "input/particle_class" in h5:
            particle_class = _decode_str_array(h5["input/particle_class"])

        lipid_mask, po4_mask, lipid_selection_mode = _select_lipid_atoms_and_po4(
            atom_names,
            molecule_ids,
            particle_class,
        )
        if not np.any(po4_mask):
            raise ValueError("No DOPC PO4 beads found in production stage file")

        lipid_indices = np.where(lipid_mask)[0]
        po4_indices = np.where(po4_mask)[0]
        n_lipid_molecules = int(np.unique(molecule_ids[po4_indices]).shape[0])
        if n_lipid_molecules < 2:
            raise ValueError("Expected at least two PO4 beads for bilayer diffusion analysis")

        reference_file = reference_stage_file if reference_stage_file is not None else stage_file
        reference_pos, reference_box_z = _load_reference_positions_and_box(reference_file)
        if reference_pos.shape[0] != n_atom:
            raise ValueError(
                f"Reference stage atom count {reference_pos.shape[0]} does not match production {n_atom}"
            )
        initial_po4_z = reference_pos[po4_indices, 2]
        top_mask, bottom_mask, z_cut_origin, z_split_threshold = _split_leaflets_from_wrapped_z(
            initial_po4_z,
            reference_box_z,
        )

        burn_start = max(1, int(math.floor(0.20 * n_frame)))
        pos = pos[burn_start:]
        times = times[burn_start:]
        box = box[burn_start:]
        n_frame = pos.shape[0]
        if n_frame < 10:
            raise ValueError("Too few frames remain after burn-in for MSD analysis")

        lipid_xy = pos[:, lipid_indices, :2]
        po4_xy = pos[:, po4_indices, :2]
        box_xy = box[:, :2]

        lipid_xy_unwrapped = _unwrap_xy(lipid_xy, box_xy)
        po4_xy_unwrapped = _unwrap_xy(po4_xy, box_xy)
        bilayer_com_xy = np.mean(lipid_xy_unwrapped, axis=1)
        drift_xy = bilayer_com_xy - bilayer_com_xy[0]
        po4_xy_corrected = po4_xy_unwrapped - drift_xy[:, None, :]

        max_lag = max(5, int(math.floor(0.25 * n_frame)))
        max_lag = min(max_lag, n_frame - 1)
        if max_lag < 5:
            raise ValueError("Not enough post-burn-in frames for lag-window analysis")

        lag_times_all, msd_all, counts_all = _compute_msd(po4_xy_corrected, times, max_lag)
        lag_times_top, msd_top, counts_top = _compute_msd(po4_xy_corrected[:, top_mask, :], times, max_lag)
        lag_times_bottom, msd_bottom, counts_bottom = _compute_msd(po4_xy_corrected[:, bottom_mask, :], times, max_lag)

        fit_all = _fit_diffusion(lag_times_all, msd_all)
        fit_top = _fit_diffusion(lag_times_top, msd_top)
        fit_bottom = _fit_diffusion(lag_times_bottom, msd_bottom)

        area_per_lipid = box[:, 0] * box[:, 1] / (0.5 * float(n_lipid_molecules))

        po4_z_wrapped = _wrapped_fraction(pos[:, po4_indices, 2], box[:, 2][:, None])
        po4_z_rotated = np.mod(po4_z_wrapped - z_cut_origin, 1.0)
        top_mean_z = np.mean(po4_z_rotated[:, top_mask], axis=1)
        bottom_mean_z = np.mean(po4_z_rotated[:, bottom_mask], axis=1)
        thickness = (top_mean_z - bottom_mean_z) * box[:, 2]

        return {
            "n_frames_total": int(burn_start + n_frame),
            "n_frames_used": int(n_frame),
            "burn_in_frames": int(burn_start),
            "lipid_selection_mode": lipid_selection_mode,
            "leaflet_assignment": "wrapped_z_gap_split",
            "leaflet_reference_file": str(reference_file),
            "leaflet_z_cut_origin_fraction": float(z_cut_origin),
            "leaflet_z_split_threshold_fraction": float(z_split_threshold),
            "n_po4": int(po4_indices.shape[0]),
            "n_lipid_molecules": int(n_lipid_molecules),
            "n_leaflet_top": int(np.sum(top_mask)),
            "n_leaflet_bottom": int(np.sum(bottom_mask)),
            "n_lipid_atoms": int(lipid_indices.shape[0]),
            "area_per_lipid_angstrom2_mean": float(np.mean(area_per_lipid)),
            "area_per_lipid_angstrom2_std": float(np.std(area_per_lipid)),
            "thickness_angstrom_mean": float(np.mean(thickness)),
            "thickness_angstrom_std": float(np.std(thickness)),
            "msd_lag_time": lag_times_all.tolist(),
            "msd_angstrom2_mean": msd_all.tolist(),
            "msd_counts": counts_all.tolist(),
            "msd_top_angstrom2_mean": msd_top.tolist(),
            "msd_top_counts": counts_top.tolist(),
            "msd_bottom_angstrom2_mean": msd_bottom.tolist(),
            "msd_bottom_counts": counts_bottom.tolist(),
            "fit_mean": fit_all,
            "fit_top": fit_top,
            "fit_bottom": fit_bottom,
            "diffusion_mean_angstrom2_per_time": float(fit_all["diffusion_angstrom2_per_time"]),
            "diffusion_mean_nm2_per_time": float(fit_all["diffusion_nm2_per_time"]),
            "diffusion_top_angstrom2_per_time": float(fit_top["diffusion_angstrom2_per_time"]),
            "diffusion_bottom_angstrom2_per_time": float(fit_bottom["diffusion_angstrom2_per_time"]),
        }


def _scale_mass_dataset(stage_file: Path, mass_scale: float) -> None:
    if abs(mass_scale - 1.0) < 1.0e-12:
        return
    with h5py.File(stage_file, "r+") as h5:
        if "input/mass" not in h5:
            raise ValueError(f"Missing input/mass in stage file: {stage_file}")
        values = np.asarray(h5["input/mass"][:], dtype=np.float64) * mass_scale
        del h5["input/mass"]
        h5.create_dataset("input/mass", data=values.astype(np.float32))


def _set_initial_position(source_file: Path, target_file: Path, env: Dict[str, str]) -> None:
    cmd = [_python_bin(), str(PREP_SCRIPT), "set-initial-position", str(source_file), str(target_file)]
    proc = sp.run(cmd, cwd=str(REPO_ROOT), env=env, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.strip() or proc.stdout.strip() or "set-initial-position failed")


def _frame_interval_time(frame_steps: int, nsteps: int, dt: float) -> float:
    effective_frame_steps = frame_steps
    if effective_frame_steps >= nsteps:
        effective_frame_steps = max(1, nsteps // 10)
    return float(effective_frame_steps) * dt


def _prepare_stage_file(
    task_dir: Path,
    manifest: Dict[str, Any],
    task: Dict[str, Any],
    stage_def: Dict[str, Any],
) -> Path:
    checkpoints_dir = task_dir / "checkpoints"
    logs_dir = task_dir / "logs"
    checkpoints_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    runtime_pdb = Path(manifest["system"]["runtime_pdb"]).expanduser().resolve()
    prepared_path = checkpoints_dir / f"{manifest['pdb_id']}.stage_{stage_def['label']}.prepared.up"

    env = _base_runtime_env()
    env["UPSIDE_NPT_ENABLE"] = "1"
    env["UPSIDE_BILAYER_LIPIDHEAD_FC"] = str(stage_def["lipidhead_fc"])

    cmd = [
        _python_bin(),
        str(PREP_SCRIPT),
        "--mode",
        "bilayer",
        "--pdb-id",
        manifest["pdb_id"],
        "--runtime-pdb-output",
        str(runtime_pdb),
        "--prepare-structure",
        "0",
        "--stage",
        stage_def["prepare_stage"],
        "--run-dir",
        str(task_dir),
        "--summary-json",
        str(task_dir / f"stage_{stage_def['label']}.summary.json"),
    ]
    _call_logged(cmd, logs_dir / f"stage_{stage_def['label']}.prepare.log", env)

    generated = task_dir / "test.input.up"
    if not generated.exists():
        raise RuntimeError(f"Stage conversion did not produce expected file: {generated}")
    if prepared_path.exists():
        prepared_path.unlink()
    generated.replace(prepared_path)
    _scale_mass_dataset(prepared_path, float(task["mass_scale"]))
    return prepared_path


def _run_minimization(stage_file: Path, task: Dict[str, Any], manifest: Dict[str, Any], log_path: Path, max_iter: int) -> None:
    dt = float(manifest["settings"]["time_step"])
    cmd = [
        str(UPSIDE_EXECUTABLE),
        str(stage_file),
        "--duration",
        "0",
        "--frame-interval",
        "1",
        "--temperature",
        str(task["temperature"]),
        "--time-step",
        str(dt),
        "--thermostat-timescale",
        str(task["damping"]),
        "--thermostat-interval",
        str(THERMOSTAT_INTERVAL),
        "--seed",
        str(task["seed"]),
        "--integrator",
        "v",
        "--disable-recentering",
        "--minimize",
        "--min-max-iter",
        str(max_iter),
        "--min-energy-tol",
        "1e-6",
        "--min-force-tol",
        "1e-3",
        "--min-step",
        "0.01",
    ]
    _call_logged(cmd, log_path, _base_runtime_env())


def _run_md(stage_file: Path, task: Dict[str, Any], manifest: Dict[str, Any], stage_def: Dict[str, Any], log_path: Path) -> None:
    dt = float(manifest["settings"]["time_step"])
    nsteps = int(stage_def["steps"])
    frame_interval = _frame_interval_time(int(stage_def["frame_steps"]), nsteps, dt)
    integrator = str(manifest["settings"].get("integrator", DEFAULT_INTEGRATOR))
    max_force = float(manifest["settings"].get("max_force", DEFAULT_MAX_FORCE))
    cmd = [
        str(UPSIDE_EXECUTABLE),
        str(stage_file),
        "--duration-steps",
        str(nsteps),
        "--frame-interval",
        str(frame_interval),
        "--temperature",
        str(task["temperature"]),
        "--time-step",
        str(dt),
        "--thermostat-timescale",
        str(task["damping"]),
        "--thermostat-interval",
        str(THERMOSTAT_INTERVAL),
        "--seed",
        str(task["seed"]),
        "--integrator",
        integrator,
        "--disable-recentering",
    ]
    if integrator == "nvtc" and max_force > 0.0:
        cmd.extend(["--max-force", str(max_force)])
    _call_logged(cmd, log_path, _base_runtime_env())


def _run_task(base_dir: Path, manifest: Dict[str, Any], task: Dict[str, Any], overwrite: bool) -> Dict[str, Any]:
    result_path = _task_result_path(base_dir, task["code"])
    if result_path.exists() and not overwrite:
        existing = _load_json(result_path)
        if existing.get("success"):
            print(f"Skipping completed task {task['code']}")
            return existing

    task_dir = _task_dir(base_dir, task["code"])
    checkpoints_dir = task_dir / "checkpoints"
    logs_dir = task_dir / "logs"
    checkpoints_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    stage_outputs: Dict[str, str] = {}
    previous_stage_file: Path | None = None

    try:
        for stage_def in _stage_schedule(manifest):
            prepared = _prepare_stage_file(task_dir, manifest, task, stage_def)
            stage_file = checkpoints_dir / f"{manifest['pdb_id']}.stage_{stage_def['label']}.up"
            shutil.copy2(prepared, stage_file)
            if previous_stage_file is not None:
                _set_initial_position(previous_stage_file, stage_file, _base_runtime_env())

            if stage_def["run_kind"] == "minimization":
                _run_minimization(
                    stage_file,
                    task,
                    manifest,
                    logs_dir / f"stage_{stage_def['label']}.run.log",
                    int(stage_def["max_iter"]),
                )
            else:
                _run_md(
                    stage_file,
                    task,
                    manifest,
                    stage_def,
                    logs_dir / f"stage_{stage_def['label']}.run.log",
                )

            stage_outputs[stage_def["label"]] = str(stage_file)
            previous_stage_file = stage_file

        production_stage = Path(stage_outputs["7.0"])
        reference_stage = Path(stage_outputs["6.0"])
        analysis = _analyze_production(production_stage, reference_stage_file=reference_stage)
        result = {
            "schema": SCHEMA_TASK_RESULT,
            "created_at_utc": _now_utc(),
            "base_dir": str(base_dir),
            "task": task,
            "success": True,
            "stage_outputs": stage_outputs,
            "analysis": analysis,
            "production_stage_file": str(production_stage),
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
            "stage_outputs": stage_outputs,
            "error": str(exc),
        }
        _write_json(result_path, result)
        raise


def _discover_stage7_analysis_tasks(base_dir: Path, manifest: Dict[str, Any]) -> List[Dict[str, Any]]:
    task_by_code = {str(task["code"]): dict(task) for task in manifest["tasks"]}
    stage7_pattern = f"{manifest['pdb_id']}.stage_7.0.up"
    discovered: List[Dict[str, Any]] = []
    for stage7_file in sorted((base_dir / "tasks").glob(f"*/checkpoints/{stage7_pattern}")):
        code = stage7_file.parents[1].name
        task = dict(task_by_code.get(code, {"code": code, "task_id": -1}))
        sort_task_id = int(task["task_id"]) if isinstance(task.get("task_id"), int) else 10**9
        reference_stage = stage7_file.with_name(f"{manifest['pdb_id']}.stage_6.0.up")
        discovered.append(
            {
                "code": code,
                "sort_task_id": sort_task_id,
                "stage7_file": str(stage7_file),
                "reference_stage_file": str(reference_stage) if reference_stage.exists() else None,
                "task": task,
            }
        )
    if not discovered:
        raise RuntimeError(f"No stage_7.0 checkpoint files found under {base_dir / 'tasks'}")

    discovered.sort(key=lambda item: (item["sort_task_id"], item["code"]))
    tasks: List[Dict[str, Any]] = []
    for analysis_task_id, item in enumerate(discovered):
        item["analysis_task_id"] = analysis_task_id
        item.pop("sort_task_id", None)
        tasks.append(item)
    return tasks


def cmd_init_stage7_analysis(args: argparse.Namespace) -> int:
    base_dir = Path(args.base_dir).expanduser().resolve()
    manifest = _load_manifest(base_dir)
    tasks = _discover_stage7_analysis_tasks(base_dir, manifest)
    payload = {
        "schema": SCHEMA_STAGE7_ANALYSIS_MANIFEST,
        "created_at_utc": _now_utc(),
        "base_dir": str(base_dir),
        "manifest_path": str(_manifest_path(base_dir)),
        "pdb_id": manifest["pdb_id"],
        "n_stage7_files": len(tasks),
        "tasks": tasks,
    }
    path = _stage7_analysis_manifest_path(base_dir)
    _write_json(path, payload)
    print(f"Stage-7 analysis manifest: {path}")
    print(f"Discovered stage_7.0 files: {len(tasks)}")
    return 0


def _run_stage7_analysis_task(
    base_dir: Path,
    analysis_manifest: Dict[str, Any],
    task_entry: Dict[str, Any],
    overwrite: bool,
) -> Dict[str, Any]:
    task = dict(task_entry["task"])
    code = str(task["code"])
    result_path = _stage7_analysis_result_path(base_dir, code)
    if result_path.exists() and not overwrite:
        existing = _load_json(result_path)
        if existing.get("success"):
            print(f"Skipping completed stage-7 analysis task {code}")
            return existing

    stage7_file = Path(task_entry["stage7_file"]).expanduser().resolve()
    if not stage7_file.exists():
        raise RuntimeError(f"Stage-7 checkpoint not found: {stage7_file}")

    reference_stage_file = None
    reference_raw = task_entry.get("reference_stage_file")
    if reference_raw:
        candidate = Path(reference_raw).expanduser().resolve()
        if candidate.exists():
            reference_stage_file = candidate

    try:
        analysis = _analyze_production(stage7_file, reference_stage_file=reference_stage_file)
        result = {
            "schema": SCHEMA_STAGE7_ANALYSIS_RESULT,
            "created_at_utc": _now_utc(),
            "base_dir": str(base_dir),
            "analysis_task_id": int(task_entry["analysis_task_id"]),
            "success": True,
            "task": task,
            "stage7_file": str(stage7_file),
            "reference_stage_file": str(reference_stage_file) if reference_stage_file is not None else None,
            "analysis_mode": "stage7_only_posthoc",
            "analysis": analysis,
        }
        _write_json(result_path, result)
        return result
    except Exception as exc:
        result = {
            "schema": SCHEMA_STAGE7_ANALYSIS_RESULT,
            "created_at_utc": _now_utc(),
            "base_dir": str(base_dir),
            "analysis_task_id": int(task_entry["analysis_task_id"]),
            "success": False,
            "task": task,
            "stage7_file": str(stage7_file),
            "reference_stage_file": str(reference_stage_file) if reference_stage_file is not None else None,
            "analysis_mode": "stage7_only_posthoc",
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


def cmd_run_stage7_analysis_local(args: argparse.Namespace) -> int:
    base_dir = Path(args.base_dir).expanduser().resolve()
    analysis_manifest = _load_stage7_analysis_manifest(base_dir)
    tasks = analysis_manifest["tasks"]
    start_task = max(0, int(args.start_task))
    selected = tasks[start_task:]
    if args.max_tasks is not None:
        selected = selected[: max(0, int(args.max_tasks))]
    for task_entry in selected:
        _run_stage7_analysis_task(base_dir, analysis_manifest, task_entry, overwrite=bool(args.overwrite))
    if not args.no_assemble:
        assemble_stage7_analysis(base_dir)
    return 0


def cmd_run_stage7_analysis_task(args: argparse.Namespace) -> int:
    round_manifest = _load_json(Path(args.round_manifest).expanduser().resolve())
    if round_manifest.get("schema") != SCHEMA_STAGE7_ANALYSIS_SLURM:
        raise RuntimeError(f"Unexpected stage-7 analysis round manifest schema: {round_manifest.get('schema')}")
    base_dir = Path(round_manifest["base_dir"]).expanduser().resolve()
    analysis_manifest = _load_stage7_analysis_manifest(base_dir)
    task_id = int(args.task_id)
    tasks = analysis_manifest["tasks"]
    if task_id < 0 or task_id >= len(tasks):
        raise RuntimeError(f"Task id out of range: {task_id} for {len(tasks)} stage-7 analysis tasks")
    _run_stage7_analysis_task(base_dir, analysis_manifest, tasks[task_id], overwrite=bool(args.overwrite))
    return 0


def _successful_results_from_dir(result_dir: Path) -> List[Dict[str, Any]]:
    results: List[Dict[str, Any]] = []
    if not result_dir.exists():
        return results
    for path in sorted(result_dir.glob("*.json")):
        payload = _load_json(path)
        if not payload.get("success"):
            continue
        if "task" not in payload or "analysis" not in payload:
            continue
        results.append(payload)
    return results


def _write_assembled_outputs(
    base_dir: Path,
    assembled_dir: Path,
    results: List[Dict[str, Any]],
    expected_replicates: int,
    total_tasks: int,
    label: str,
) -> int:
    assembled_dir.mkdir(parents=True, exist_ok=True)

    task_rows: List[Dict[str, Any]] = []
    grouped: Dict[Tuple[float, float, float], List[Dict[str, Any]]] = defaultdict(list)
    for result in results:
        task = result["task"]
        analysis = result["analysis"]
        row = {
            "task_id": int(task["task_id"]),
            "code": task["code"],
            "temperature": float(task["temperature"]),
            "damping": float(task["damping"]),
            "mass_scale": float(task["mass_scale"]),
            "replicate": int(task["replicate"]),
            "diffusion_mean_angstrom2_per_time": float(analysis["diffusion_mean_angstrom2_per_time"]),
            "diffusion_mean_nm2_per_time": float(analysis["diffusion_mean_nm2_per_time"]),
            "diffusion_top_angstrom2_per_time": float(analysis["diffusion_top_angstrom2_per_time"]),
            "diffusion_bottom_angstrom2_per_time": float(analysis["diffusion_bottom_angstrom2_per_time"]),
            "area_per_lipid_angstrom2_mean": float(analysis["area_per_lipid_angstrom2_mean"]),
            "thickness_angstrom_mean": float(analysis["thickness_angstrom_mean"]),
            "fit_r2_mean": float(analysis["fit_mean"]["r2"]),
            "fit_r2_top": float(analysis["fit_top"]["r2"]),
            "fit_r2_bottom": float(analysis["fit_bottom"]["r2"]),
        }
        task_rows.append(row)
        grouped[(row["temperature"], row["damping"], row["mass_scale"])].append(row)

    task_csv = assembled_dir / "task_results.csv"
    with task_csv.open("w", encoding="utf-8", newline="") as fh:
        fieldnames = [
            "task_id",
            "code",
            "temperature",
            "damping",
            "mass_scale",
            "replicate",
            "diffusion_mean_angstrom2_per_time",
            "diffusion_mean_nm2_per_time",
            "diffusion_top_angstrom2_per_time",
            "diffusion_bottom_angstrom2_per_time",
            "area_per_lipid_angstrom2_mean",
            "thickness_angstrom_mean",
            "fit_r2_mean",
            "fit_r2_top",
            "fit_r2_bottom",
        ]
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in sorted(task_rows, key=lambda item: item["task_id"]):
            writer.writerow(row)

    summary_rows: List[Dict[str, Any]] = []
    for (temperature, damping, mass_scale), rows in sorted(grouped.items()):
        d_values = np.asarray([row["diffusion_mean_angstrom2_per_time"] for row in rows], dtype=np.float64)
        area_values = np.asarray([row["area_per_lipid_angstrom2_mean"] for row in rows], dtype=np.float64)
        thickness_values = np.asarray([row["thickness_angstrom_mean"] for row in rows], dtype=np.float64)
        summary_rows.append(
            {
                "temperature": float(temperature),
                "damping": float(damping),
                "mass_scale": float(mass_scale),
                "n_replicates_expected": int(expected_replicates),
                "n_replicates_completed": int(len(rows)),
                "diffusion_mean_angstrom2_per_time_mean": float(np.mean(d_values)),
                "diffusion_mean_angstrom2_per_time_std": float(np.std(d_values)),
                "diffusion_mean_nm2_per_time_mean": float(np.mean(d_values) / 100.0),
                "diffusion_mean_nm2_per_time_std": float(np.std(d_values) / 100.0),
                "area_per_lipid_angstrom2_mean": float(np.mean(area_values)),
                "area_per_lipid_angstrom2_std": float(np.std(area_values)),
                "thickness_angstrom_mean": float(np.mean(thickness_values)),
                "thickness_angstrom_std": float(np.std(thickness_values)),
            }
        )

    summary_csv = assembled_dir / "condition_summary.csv"
    with summary_csv.open("w", encoding="utf-8", newline="") as fh:
        fieldnames = [
            "temperature",
            "damping",
            "mass_scale",
            "n_replicates_expected",
            "n_replicates_completed",
            "diffusion_mean_angstrom2_per_time_mean",
            "diffusion_mean_angstrom2_per_time_std",
            "diffusion_mean_nm2_per_time_mean",
            "diffusion_mean_nm2_per_time_std",
            "area_per_lipid_angstrom2_mean",
            "area_per_lipid_angstrom2_std",
            "thickness_angstrom_mean",
            "thickness_angstrom_std",
        ]
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in summary_rows:
            writer.writerow(row)

    payload = {
        "created_at_utc": _now_utc(),
        "base_dir": str(base_dir),
        "result_label": label,
        "n_tasks_total": int(total_tasks),
        "n_tasks_completed_successfully": len(task_rows),
        "n_conditions_completed": len(summary_rows),
        "task_results_csv": str(task_csv),
        "condition_summary_csv": str(summary_csv),
    }
    _write_json(assembled_dir / "summary.json", payload)
    print(f"Assembled {label} results written under {assembled_dir}")
    return 0


def _successful_results(base_dir: Path) -> List[Dict[str, Any]]:
    return _successful_results_from_dir(_task_results_dir(base_dir))


def assemble_results(base_dir: Path) -> int:
    manifest = _load_manifest(base_dir)
    return _write_assembled_outputs(
        base_dir,
        _assembled_dir(base_dir),
        _successful_results(base_dir),
        int(manifest["settings"]["replicates"]),
        len(manifest["tasks"]),
        "simulation",
    )


def assemble_stage7_analysis(base_dir: Path) -> int:
    manifest = _load_manifest(base_dir)
    analysis_manifest = _load_stage7_analysis_manifest(base_dir)
    return _write_assembled_outputs(
        base_dir,
        _stage7_analysis_assembled_dir(base_dir),
        _successful_results_from_dir(_stage7_analysis_results_dir(base_dir)),
        int(manifest["settings"]["replicates"]),
        len(analysis_manifest["tasks"]),
        "stage7-analysis",
    )


def cmd_assemble_stage7_analysis(args: argparse.Namespace) -> int:
    base_dir = Path(args.base_dir).expanduser().resolve()
    return assemble_stage7_analysis(base_dir)


def _stage7_analysis_array_script_content(round_manifest_path: Path, base_dir: Path, n_tasks: int) -> str:
    output_path = _stage7_analysis_slurm_dir(base_dir) / "analyze-%A_%a.out"
    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name=bilayer-diff-a7-{base_dir.name}",
        f"#SBATCH --output={output_path}",
        f"#SBATCH --time={os.environ.get('BILAYER_DIFF_ANALYSIS_WALLTIME', '08:00:00')}",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        f"#SBATCH --cpus-per-task={int(os.environ.get('BILAYER_DIFF_ANALYSIS_CPUS_PER_TASK', '1'))}",
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
            '  module load "${BILAYER_DIFF_HDF5_MODULE:-hdf5/1.14.3}" || true',
            "fi",
            "",
            'if [ -f "$SCRIPT_DIR/.venv/bin/activate" ]; then',
            '  source "$SCRIPT_DIR/.venv/bin/activate"',
            'elif [ -f "$PROJECT_ROOT/.venv/bin/activate" ]; then',
            '  source "$PROJECT_ROOT/.venv/bin/activate"',
            "fi",
            "",
            'if [ -n "${BILAYER_DIFF_PYTHON:-}" ]; then',
            '  if command -v "${BILAYER_DIFF_PYTHON}" >/dev/null 2>&1; then',
            '    PYTHON_BIN="$(command -v "${BILAYER_DIFF_PYTHON}")"',
            "  else",
            '    PYTHON_BIN="${BILAYER_DIFF_PYTHON}"',
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
            'cmd=("$PYTHON_BIN" -u "$SCRIPT_DIR/workflow.py" run-stage7-analysis-task --round-manifest "$ROUND_MANIFEST" --task-id "$SLURM_ARRAY_TASK_ID")',
            'if [ "${BILAYER_DIFF_OVERWRITE:-0}" = "1" ]; then',
            '  cmd+=(--overwrite)',
            "fi",
            '"${cmd[@]}"',
        ]
    )
    return "\n".join(lines) + "\n"


def _stage7_analysis_collector_script_content(round_manifest_path: Path, base_dir: Path) -> str:
    output_path = _stage7_analysis_slurm_dir(base_dir) / "collect-%j.out"
    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name=bilayer-diff-a7-collect-{base_dir.name}",
        f"#SBATCH --output={output_path}",
        f"#SBATCH --time={os.environ.get('BILAYER_DIFF_ANALYSIS_COLLECT_WALLTIME', '01:00:00')}",
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
            '  module load "${BILAYER_DIFF_HDF5_MODULE:-hdf5/1.14.3}" || true',
            "fi",
            "",
            'if [ -f "$SCRIPT_DIR/.venv/bin/activate" ]; then',
            '  source "$SCRIPT_DIR/.venv/bin/activate"',
            'elif [ -f "$PROJECT_ROOT/.venv/bin/activate" ]; then',
            '  source "$PROJECT_ROOT/.venv/bin/activate"',
            "fi",
            'if [ -n "${BILAYER_DIFF_PYTHON:-}" ]; then',
            '  if command -v "${BILAYER_DIFF_PYTHON}" >/dev/null 2>&1; then',
            '    PYTHON_BIN="$(command -v "${BILAYER_DIFF_PYTHON}")"',
            "  else",
            '    PYTHON_BIN="${BILAYER_DIFF_PYTHON}"',
            "  fi",
            'elif [ -n "${VIRTUAL_ENV:-}" ] && [ -x "${VIRTUAL_ENV}/bin/python3" ]; then',
            '  PYTHON_BIN="${VIRTUAL_ENV}/bin/python3"',
            "else",
            '  PYTHON_BIN="$(command -v python3)"',
            "fi",
            'export PYTHONUNBUFFERED=1',
            '"$PYTHON_BIN" -u "$SCRIPT_DIR/workflow.py" assemble-stage7-analysis --base-dir "$BASE_DIR"',
        ]
    )
    return "\n".join(lines) + "\n"


def cmd_submit_stage7_analysis_slurm(args: argparse.Namespace) -> int:
    base_dir = Path(args.base_dir).expanduser().resolve()
    analysis_manifest = _load_stage7_analysis_manifest(base_dir)
    n_tasks = len(analysis_manifest["tasks"])
    if n_tasks <= 0:
        raise RuntimeError("Stage-7 analysis manifest contains no tasks")

    round_manifest = {
        "schema": SCHEMA_STAGE7_ANALYSIS_SLURM,
        "created_at_utc": _now_utc(),
        "base_dir": str(base_dir),
        "analysis_manifest_path": str(_stage7_analysis_manifest_path(base_dir)),
        "n_tasks": n_tasks,
    }
    round_manifest_path = _stage7_analysis_slurm_dir(base_dir) / "round_manifest.json"
    _write_json(round_manifest_path, round_manifest)

    array_script = _stage7_analysis_slurm_dir(base_dir) / "analyze_stage7_array.sbatch"
    collector_script = _stage7_analysis_slurm_dir(base_dir) / "collect_stage7_analysis.sbatch"
    _write_text(
        array_script,
        _stage7_analysis_array_script_content(round_manifest_path, base_dir, n_tasks),
        executable=True,
    )
    _write_text(
        collector_script,
        _stage7_analysis_collector_script_content(round_manifest_path, base_dir),
        executable=True,
    )

    print(f"Stage-7 analysis manifest: {_stage7_analysis_manifest_path(base_dir)}")
    print(f"Stage-7 round manifest: {round_manifest_path}")
    print(f"Stage-7 array script: {array_script}")
    print(f"Stage-7 collector script: {collector_script}")

    if args.no_submit:
        print("Staged stage-7 Slurm scripts only; no jobs submitted.")
        return 0

    if not _shutil_which("sbatch"):
        raise RuntimeError("sbatch not available in PATH")

    analyze_job = _run_sbatch(["sbatch", "--parsable", str(array_script)])
    collect_job = _run_sbatch(["sbatch", "--parsable", f"--dependency=afterok:{analyze_job}", str(collector_script)])
    _write_json(
        _stage7_analysis_slurm_dir(base_dir) / "submission.json",
        {
            "created_at_utc": _now_utc(),
            "analyze_job_id": analyze_job,
            "collect_job_id": collect_job,
        },
    )
    print(f"Submitted stage-7 analysis array job: {analyze_job}")
    print(f"Submitted stage-7 analysis collector job: {collect_job}")
    return 0


def cmd_pack_stage7_results(args: argparse.Namespace) -> int:
    base_dir = Path(args.base_dir).expanduser().resolve()
    output = (
        Path(args.output).expanduser().resolve()
        if args.output
        else _download_dir(base_dir) / f"{base_dir.name}-stage7-results.tar.gz"
    )
    include_stage7_files = not bool(args.metadata_only)

    archive_sources: List[Path] = []
    for path in [
        _manifest_path(base_dir),
        _assembled_dir(base_dir),
        _results_dir(base_dir),
        _stage7_analysis_dir(base_dir),
    ]:
        if path.exists():
            archive_sources.append(path)

    stage7_files: List[Path] = []
    if include_stage7_files:
        stage7_files = sorted(base_dir.glob(f"tasks/*/checkpoints/{BASE_PDB_ID}.stage_7.0.up"))
        archive_sources.extend(stage7_files)

    if not archive_sources:
        raise RuntimeError(f"No files available to pack under {base_dir}")

    output.parent.mkdir(parents=True, exist_ok=True)
    archive_root = Path(base_dir.name)
    with tarfile.open(output, "w:gz") as tar:
        for path in archive_sources:
            tar.add(path, arcname=str(archive_root / path.relative_to(base_dir)), recursive=True)

    manifest_path = Path(str(output) + ".manifest.json")
    _write_json(
        manifest_path,
        {
            "created_at_utc": _now_utc(),
            "base_dir": str(base_dir),
            "output_archive": str(output),
            "output_manifest": str(manifest_path),
            "include_stage7_files": include_stage7_files,
            "n_archive_sources": len(archive_sources),
            "n_stage7_files": len(stage7_files),
        },
    )
    print(f"Created archive: {output}")
    print(f"Archive manifest: {manifest_path}")
    return 0


def cmd_assemble_results(args: argparse.Namespace) -> int:
    base_dir = Path(args.base_dir).expanduser().resolve()
    return assemble_results(base_dir)


def _array_script_content(round_manifest_path: Path, base_dir: Path, n_tasks: int) -> str:
    output_path = _slurm_dir(base_dir) / "train-%A_%a.out"
    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name=bilayer-diff-{base_dir.name}",
        f"#SBATCH --output={output_path}",
        f"#SBATCH --time={os.environ.get('BILAYER_DIFF_TRAIN_WALLTIME', '24:00:00')}",
        "#SBATCH --nodes=1",
        "#SBATCH --ntasks=1",
        f"#SBATCH --cpus-per-task={int(os.environ.get('BILAYER_DIFF_CPUS_PER_TASK', '1'))}",
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
            '  module load "${BILAYER_DIFF_HDF5_MODULE:-hdf5/1.14.3}" || true',
            "fi",
            "",
            'if [ -f "$SCRIPT_DIR/.venv/bin/activate" ]; then',
            '  source "$SCRIPT_DIR/.venv/bin/activate"',
            'elif [ -f "$PROJECT_ROOT/.venv/bin/activate" ]; then',
            '  source "$PROJECT_ROOT/.venv/bin/activate"',
            "fi",
            "",
            'if [ -n "${BILAYER_DIFF_PYTHON:-}" ]; then',
            '  if command -v "${BILAYER_DIFF_PYTHON}" >/dev/null 2>&1; then',
            '    PYTHON_BIN="$(command -v "${BILAYER_DIFF_PYTHON}")"',
            "  else",
            '    PYTHON_BIN="${BILAYER_DIFF_PYTHON}"',
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
            '"$PYTHON_BIN" -u "$SCRIPT_DIR/workflow.py" run-array-task --round-manifest "$ROUND_MANIFEST" --task-id "$SLURM_ARRAY_TASK_ID"',
        ]
    )
    return "\n".join(lines) + "\n"


def _collector_script_content(round_manifest_path: Path, base_dir: Path) -> str:
    output_path = _slurm_dir(base_dir) / "collect-%j.out"
    lines = [
        "#!/bin/bash",
        f"#SBATCH --job-name=bilayer-diff-collect-{base_dir.name}",
        f"#SBATCH --output={output_path}",
        f"#SBATCH --time={os.environ.get('BILAYER_DIFF_COLLECT_WALLTIME', '01:00:00')}",
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
            '  module load "${BILAYER_DIFF_HDF5_MODULE:-hdf5/1.14.3}" || true',
            "fi",
            "",
            'if [ -f "$SCRIPT_DIR/.venv/bin/activate" ]; then',
            '  source "$SCRIPT_DIR/.venv/bin/activate"',
            'elif [ -f "$PROJECT_ROOT/.venv/bin/activate" ]; then',
            '  source "$PROJECT_ROOT/.venv/bin/activate"',
            "fi",
            'if [ -n "${BILAYER_DIFF_PYTHON:-}" ]; then',
            '  if command -v "${BILAYER_DIFF_PYTHON}" >/dev/null 2>&1; then',
            '    PYTHON_BIN="$(command -v "${BILAYER_DIFF_PYTHON}")"',
            "  else",
            '    PYTHON_BIN="${BILAYER_DIFF_PYTHON}"',
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
    _write_text(collector_script, _collector_script_content(round_manifest_path, base_dir), executable=True)

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
        _slurm_dir(base_dir) / "round_submission.json",
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
    parser = argparse.ArgumentParser(description="Bilayer lateral diffusion workflow")
    sub = parser.add_subparsers(dest="command", required=True)

    p_init = sub.add_parser("init-run", help="Initialize a new bilayer diffusion run")
    p_init.add_argument("--base-dir", default=str(WORKFLOW_DIR / "runs" / "default"))
    p_init.add_argument("--bilayer-pdb", default=str(DEFAULT_BILAYER_PDB))
    p_init.add_argument("--damping-values", type=_float_list_arg, default=DEFAULT_DAMPING_VALUES)
    p_init.add_argument("--mass-scales", type=_float_list_arg, default=DEFAULT_MASS_SCALES)
    p_init.add_argument("--temperature-values", type=_float_list_arg, default=DEFAULT_TEMPERATURE_VALUES)
    p_init.add_argument("--replicates", type=int, default=DEFAULT_REPLICATES)
    p_init.add_argument("--seed", type=int, default=DEFAULT_SEED)
    p_init.add_argument("--xy-scale", type=float, default=DEFAULT_XY_SCALE)
    p_init.add_argument("--box-padding-z", type=float, default=DEFAULT_BOX_PADDING_Z)
    p_init.add_argument("--salt-molar", type=float, default=DEFAULT_SALT_MOLAR)
    p_init.add_argument("--ion-cutoff", type=float, default=DEFAULT_ION_CUTOFF)
    p_init.add_argument("--min-60-max-iter", type=int, default=DEFAULT_MIN_60_MAX_ITER)
    p_init.add_argument("--min-61-max-iter", type=int, default=DEFAULT_MIN_61_MAX_ITER)
    p_init.add_argument("--equilibration-steps", type=int, default=DEFAULT_EQ_STEPS)
    p_init.add_argument("--production-steps", type=int, default=DEFAULT_PROD_STEPS)
    p_init.add_argument("--time-step", type=float, default=DEFAULT_TIME_STEP)
    p_init.add_argument("--integrator", choices=["v", "mv", "nvtc"], default=DEFAULT_INTEGRATOR)
    p_init.add_argument("--max-force", type=float, default=DEFAULT_MAX_FORCE)
    p_init.add_argument("--equilibration-frame-steps", type=int, default=DEFAULT_EQ_FRAME_STEPS)
    p_init.add_argument("--production-frame-steps", type=int, default=DEFAULT_PROD_FRAME_STEPS)
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

    p_init_stage7 = sub.add_parser("init-stage7-analysis", help="Discover existing stage-7 checkpoint files for post-analysis")
    p_init_stage7.add_argument("--base-dir", default=str(WORKFLOW_DIR / "runs" / "default"))
    p_init_stage7.set_defaults(func=cmd_init_stage7_analysis)

    p_stage7_local = sub.add_parser("run-stage7-analysis-local", help="Run stage-7 post-analysis tasks locally")
    p_stage7_local.add_argument("--base-dir", default=str(WORKFLOW_DIR / "runs" / "default"))
    p_stage7_local.add_argument("--max-tasks", type=int, default=None)
    p_stage7_local.add_argument("--start-task", type=int, default=0)
    p_stage7_local.add_argument("--overwrite", action="store_true")
    p_stage7_local.add_argument("--no-assemble", action="store_true")
    p_stage7_local.set_defaults(func=cmd_run_stage7_analysis_local)

    p_stage7_task = sub.add_parser("run-stage7-analysis-task", help="Run one stage-7 post-analysis task")
    p_stage7_task.add_argument("--round-manifest", required=True)
    p_stage7_task.add_argument("--task-id", required=True, type=int)
    p_stage7_task.add_argument("--overwrite", action="store_true")
    p_stage7_task.set_defaults(func=cmd_run_stage7_analysis_task)

    p_stage7_assemble = sub.add_parser("assemble-stage7-analysis", help="Aggregate completed stage-7 post-analysis results")
    p_stage7_assemble.add_argument("--base-dir", default=str(WORKFLOW_DIR / "runs" / "default"))
    p_stage7_assemble.set_defaults(func=cmd_assemble_stage7_analysis)

    p_stage7_submit = sub.add_parser("submit-stage7-analysis-slurm", help="Stage and optionally submit stage-7 analysis Slurm scripts")
    p_stage7_submit.add_argument("--base-dir", default=str(WORKFLOW_DIR / "runs" / "default"))
    p_stage7_submit.add_argument("--no-submit", action="store_true")
    p_stage7_submit.set_defaults(func=cmd_submit_stage7_analysis_slurm)

    p_pack = sub.add_parser("pack-stage7-results", help="Pack stage-7 results and analysis outputs for download")
    p_pack.add_argument("--base-dir", default=str(WORKFLOW_DIR / "runs" / "default"))
    p_pack.add_argument("--output", default=None)
    p_pack.add_argument("--metadata-only", action="store_true")
    p_pack.set_defaults(func=cmd_pack_stage7_results)

    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
