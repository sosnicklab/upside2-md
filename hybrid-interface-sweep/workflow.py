#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import shlex
import subprocess as sp
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Sequence

import h5py
import numpy as np


WORKFLOW_DIR = Path(os.path.abspath(__file__)).parent
REPO_ROOT = WORKFLOW_DIR.parent
HYBRID_RUN_SCRIPT = REPO_ROOT / "example" / "16.MARTINI" / "run_sim_1rkl.sh"

SCHEMA_MANIFEST = "hybrid_interface_sweep_manifest_v1"
SCHEMA_SLURM_ROUND = "hybrid_interface_sweep_slurm_round_v1"
SCHEMA_TASK_RESULT = "hybrid_interface_sweep_task_result_v1"
SCHEMA_ANALYSIS_MANIFEST = "hybrid_interface_sweep_analysis_manifest_v1"
SCHEMA_ANALYSIS_SLURM = "hybrid_interface_sweep_analysis_slurm_v1"
SCHEMA_ANALYSIS_RESULT = "hybrid_interface_sweep_analysis_result_v1"

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
    molecule_ids: np.ndarray,
    particle_class: np.ndarray | None,
    protein_membership: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    po4_mask = (atom_names == "PO4") & (protein_membership < 0)
    if particle_class is not None:
        lipid_mask = particle_class == "LIPID"
        lipid_po4_mask = lipid_mask & po4_mask
        if np.any(lipid_po4_mask):
            return lipid_mask, lipid_po4_mask

    if not np.any(po4_mask):
        raise ValueError("No non-protein PO4 beads found in stage file")

    lipid_molecule_ids = np.unique(molecule_ids[po4_mask])
    lipid_mask = np.isin(molecule_ids, lipid_molecule_ids) & (protein_membership < 0)
    lipid_po4_mask = lipid_mask & (atom_names == "PO4")
    if not np.any(lipid_po4_mask):
        raise ValueError("Could not infer lipid molecules from PO4 beads")
    return lipid_mask, lipid_po4_mask


def _select_protein_ca_indices(h5: h5py.File, atom_names: np.ndarray, protein_membership: np.ndarray) -> np.ndarray:
    inp = h5["input"]
    if "hybrid_bb_map" in inp:
        bb = inp["hybrid_bb_map"]
        atom_indices = np.asarray(bb["atom_indices"][:], dtype=np.int32)
        atom_mask = np.asarray(bb["atom_mask"][:], dtype=np.int32)
        if atom_indices.ndim == 2 and atom_indices.shape[1] >= 2:
            ca_mask = atom_mask[:, 1] != 0
            ca_indices = atom_indices[ca_mask, 1]
            ca_indices = ca_indices[ca_indices >= 0]
            if ca_indices.size:
                return np.unique(ca_indices.astype(np.int32, copy=False))

    fallback = np.where((protein_membership >= 0) & (atom_names == "CA"))[0].astype(np.int32, copy=False)
    if fallback.size == 0:
        raise ValueError("Could not identify protein CA carrier atoms for analysis")
    return fallback


def _analyze_stage7(stage_file: Path) -> Dict[str, Any]:
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
        molecule_ids = np.asarray(inp["molecule_ids"][:], dtype=np.int32)
        particle_class = _decode_str_array(inp["particle_class"]) if "particle_class" in inp else None
        protein_membership = np.asarray(inp["hybrid_env_topology"]["protein_membership"][:], dtype=np.int32)
        if protein_membership.shape[0] != n_atom:
            raise ValueError("protein_membership length does not match atom count")

        lipid_mask, po4_mask = _select_lipid_atoms_and_po4(atom_names, molecule_ids, particle_class, protein_membership)
        protein_ca_indices = _select_protein_ca_indices(h5, atom_names, protein_membership)
        po4_indices = np.where(po4_mask)[0]
        lipid_indices = np.where(lipid_mask)[0]
        if lipid_indices.size == 0:
            raise ValueError("No lipid atoms selected for bilayer COM analysis")

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
        protein_ca_xy = pos[:, protein_ca_indices, :2]

        lipid_xy_unwrapped = _unwrap_xy(lipid_xy, box_xy)
        po4_xy_unwrapped = _unwrap_xy(po4_xy, box_xy)
        protein_ca_xy_unwrapped = _unwrap_xy(protein_ca_xy, box_xy)

        bilayer_com_xy = np.mean(lipid_xy_unwrapped, axis=1)
        po4_xy_corrected = po4_xy_unwrapped - bilayer_com_xy[:, None, :]
        protein_com_xy = np.mean(protein_ca_xy_unwrapped, axis=1)
        protein_rel_xy = protein_com_xy - bilayer_com_xy

        max_lag = max(5, int(math.floor(0.25 * n_frame)))
        max_lag = min(max_lag, n_frame - 1)
        if max_lag < 5:
            raise ValueError("Not enough post-burn-in frames for lag-window analysis")

        lag_times_protein, msd_protein, counts_protein = _compute_msd(
            protein_rel_xy[:, None, :], times, max_lag
        )
        lag_times_po4, msd_po4, counts_po4 = _compute_msd(po4_xy_corrected, times, max_lag)
        fit_protein = _fit_diffusion(lag_times_protein, msd_protein)
        fit_po4 = _fit_diffusion(lag_times_po4, msd_po4)

        return {
            "n_frames_total": int(burn_start + n_frame),
            "n_frames_used": int(n_frame),
            "burn_in_frames": int(burn_start),
            "n_atom": int(n_atom),
            "n_lipid_atoms": int(lipid_indices.shape[0]),
            "n_po4": int(po4_indices.shape[0]),
            "n_protein_ca": int(protein_ca_indices.shape[0]),
            "protein_ca_indices": protein_ca_indices.astype(int).tolist(),
            "msd_lag_time": lag_times_protein.tolist(),
            "protein_com_msd_angstrom2": msd_protein.tolist(),
            "protein_com_msd_counts": counts_protein.tolist(),
            "po4_msd_angstrom2": msd_po4.tolist(),
            "po4_msd_counts": counts_po4.tolist(),
            "fit_protein": fit_protein,
            "fit_po4": fit_po4,
            "protein_lateral_diffusion_angstrom2_per_time": float(fit_protein["diffusion_angstrom2_per_time"]),
            "protein_lateral_diffusion_nm2_per_time": float(fit_protein["diffusion_nm2_per_time"]),
            "po4_lateral_diffusion_angstrom2_per_time": float(fit_po4["diffusion_angstrom2_per_time"]),
            "po4_lateral_diffusion_nm2_per_time": float(fit_po4["diffusion_nm2_per_time"]),
        }


def _discover_analysis_tasks(base_dir: Path, manifest: Dict[str, Any]) -> List[Dict[str, Any]]:
    task_by_code = {str(task["code"]): dict(task) for task in manifest["tasks"]}
    discovered: List[Dict[str, Any]] = []
    pattern = f"{manifest['settings']['pdb_id']}.stage_7.0.up"
    for stage7_file in sorted(base_dir.glob(f"tasks/*/run/checkpoints/{pattern}")):
        code = stage7_file.parents[2].name
        task = dict(task_by_code.get(code, {"code": code, "task_id": -1}))
        sort_task_id = int(task["task_id"]) if isinstance(task.get("task_id"), int) else 10**9
        discovered.append(
            {
                "code": code,
                "sort_task_id": sort_task_id,
                "stage7_file": str(stage7_file),
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


def _run_analysis_task(
    base_dir: Path,
    task_entry: Dict[str, Any],
    overwrite: bool,
) -> Dict[str, Any]:
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
        analysis = _analyze_stage7(stage7_file)
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


def assemble_analysis(base_dir: Path) -> int:
    manifest = _load_manifest(base_dir)
    analysis_manifest = _load_analysis_manifest(base_dir)
    results_dir = _analysis_results_dir(base_dir)
    assembled_dir = _analysis_assembled_dir(base_dir)
    assembled_dir.mkdir(parents=True, exist_ok=True)

    successful_results = _successful_results_from_dir(results_dir)
    task_rows: List[Dict[str, Any]] = []
    grouped: Dict[float, List[Dict[str, Any]]] = {}
    for result in successful_results:
        task = result["task"]
        analysis = result["analysis"]
        row = {
            "analysis_task_id": int(result["analysis_task_id"]),
            "task_id": int(task.get("task_id", -1)),
            "code": str(task["code"]),
            "interface_scale": float(task["interface_scale"]),
            "replicate": int(task["replicate"]),
            "protein_lateral_diffusion_angstrom2_per_time": float(
                analysis["protein_lateral_diffusion_angstrom2_per_time"]
            ),
            "protein_lateral_diffusion_nm2_per_time": float(
                analysis["protein_lateral_diffusion_nm2_per_time"]
            ),
            "po4_lateral_diffusion_angstrom2_per_time": float(
                analysis["po4_lateral_diffusion_angstrom2_per_time"]
            ),
            "po4_lateral_diffusion_nm2_per_time": float(
                analysis["po4_lateral_diffusion_nm2_per_time"]
            ),
            "protein_fit_r2": float(analysis["fit_protein"]["r2"]),
            "po4_fit_r2": float(analysis["fit_po4"]["r2"]),
            "n_frames_used": int(analysis["n_frames_used"]),
            "n_protein_ca": int(analysis["n_protein_ca"]),
            "n_po4": int(analysis["n_po4"]),
            "stage7_file": str(result["stage7_file"]),
        }
        task_rows.append(row)
        grouped.setdefault(row["interface_scale"], []).append(row)

    task_csv = assembled_dir / "task_results.csv"
    with task_csv.open("w", encoding="utf-8", newline="") as fh:
        fieldnames = [
            "analysis_task_id",
            "task_id",
            "code",
            "interface_scale",
            "replicate",
            "protein_lateral_diffusion_angstrom2_per_time",
            "protein_lateral_diffusion_nm2_per_time",
            "po4_lateral_diffusion_angstrom2_per_time",
            "po4_lateral_diffusion_nm2_per_time",
            "protein_fit_r2",
            "po4_fit_r2",
            "n_frames_used",
            "n_protein_ca",
            "n_po4",
            "stage7_file",
        ]
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in sorted(task_rows, key=lambda item: (item["analysis_task_id"], item["task_id"], item["code"])):
            writer.writerow(row)

    expected_replicates = int(manifest["settings"]["replicates"])
    summary_rows: List[Dict[str, Any]] = []
    for interface_scale, rows in sorted(grouped.items(), reverse=True):
        protein_nm = np.asarray(
            [row["protein_lateral_diffusion_nm2_per_time"] for row in rows],
            dtype=np.float64,
        )
        po4_nm = np.asarray(
            [row["po4_lateral_diffusion_nm2_per_time"] for row in rows],
            dtype=np.float64,
        )
        protein_r2 = np.asarray([row["protein_fit_r2"] for row in rows], dtype=np.float64)
        po4_r2 = np.asarray([row["po4_fit_r2"] for row in rows], dtype=np.float64)
        summary_rows.append(
            {
                "interface_scale": float(interface_scale),
                "n_replicates_expected": expected_replicates,
                "n_replicates_completed": len(rows),
                "protein_diffusion_nm2_per_time_mean": float(np.mean(protein_nm)),
                "protein_diffusion_nm2_per_time_std": float(np.std(protein_nm, ddof=0)),
                "po4_diffusion_nm2_per_time_mean": float(np.mean(po4_nm)),
                "po4_diffusion_nm2_per_time_std": float(np.std(po4_nm, ddof=0)),
                "protein_fit_r2_min": float(np.min(protein_r2)),
                "po4_fit_r2_min": float(np.min(po4_r2)),
            }
        )

    summary_csv = assembled_dir / "condition_summary.csv"
    with summary_csv.open("w", encoding="utf-8", newline="") as fh:
        fieldnames = [
            "interface_scale",
            "n_replicates_expected",
            "n_replicates_completed",
            "protein_diffusion_nm2_per_time_mean",
            "protein_diffusion_nm2_per_time_std",
            "po4_diffusion_nm2_per_time_mean",
            "po4_diffusion_nm2_per_time_std",
            "protein_fit_r2_min",
            "po4_fit_r2_min",
        ]
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
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

    payload = {
        "created_at_utc": _now_utc(),
        "base_dir": str(base_dir),
        "n_analysis_tasks_total": int(len(analysis_manifest["tasks"])),
        "n_analysis_tasks_completed_successfully": len(task_rows),
        "n_analysis_tasks_failed": len(failed_results),
        "task_results_csv": str(task_csv),
        "condition_summary_csv": str(summary_csv),
        "failed_tasks_csv": str(failed_csv),
    }
    _write_json(assembled_dir / "summary.json", payload)
    print(f"Assembled analysis results written under {assembled_dir}")
    return 0


def cmd_assemble_analysis(args: argparse.Namespace) -> int:
    base_dir = Path(args.base_dir).expanduser().resolve()
    return assemble_analysis(base_dir)


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

    array_script = _slurm_dir(base_dir) / "train_array.sbatch"
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
            "train_job_id": train_job,
            "collect_job_id": collect_job,
        },
    )
    print(f"Submitted training array job: {train_job}")
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
