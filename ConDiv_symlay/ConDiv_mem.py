#!/usr/bin/env python3
"""Modernized membrane ConDiv workflow (Python 3, local + Slurm compatible)."""

from __future__ import annotations

import collections
import itertools
import json
import multiprocessing as mp
import os
import pickle as cp
import re
import shutil
import socket
import subprocess as sp
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import mdtraj as md
import numpy as np
import tables as tb


# Keep module name stable for pickling even when executed as a script.
sys.modules.setdefault("ConDiv_mem", sys.modules[__name__])


def _locate_project_root() -> Path:
    env_root = os.environ.get("CONDIV_PROJECT_ROOT")
    if env_root:
        root = Path(env_root).expanduser().resolve()
        if (root / "py" / "run_upside.py").exists():
            return root

    starts = [Path(__file__).resolve(), Path.cwd().resolve()]
    seen: set[Path] = set()
    for start in starts:
        for candidate in [start.parent, *start.parents]:
            if candidate in seen:
                continue
            seen.add(candidate)
            if (candidate / "py" / "run_upside.py").exists() and (candidate / "parameters").exists():
                return candidate

    raise RuntimeError(
        "Could not locate project root containing py/run_upside.py. "
        "Set CONDIV_PROJECT_ROOT to /Users/yinhan/Documents/upside2-md."
    )


PROJECT_ROOT = _locate_project_root()
PY_DIR = PROJECT_ROOT / "py"
OBJ_DIR = PROJECT_ROOT / "obj"
if str(PY_DIR) not in sys.path:
    sys.path.insert(0, str(PY_DIR))
if str(OBJ_DIR) not in sys.path:
    sys.path.insert(0, str(OBJ_DIR))

import mdtraj_upside as mu
import run_upside as ru
import upside_engine as ue
from symlay_utils import (
    PAIR_TEACHER_SCHEMA,
    DEFAULT_DENSE_GRID_SIZE,
    DEFAULT_SUPPORT_MARGIN,
    build_pair_teacher,
    compute_pair_regularization_gradient,
    load_layer_manifest,
    load_membrane_arrays,
    project_membrane_arrays,
    repo_relative,
    write_projected_membrane_file,
)

np.set_printoptions(precision=3, suppress=True)


def _env_bool(name: str, default: bool) -> bool:
    raw = os.environ.get(name)
    if raw is None:
        return default
    return raw.strip().lower() in {"1", "true", "yes", "on"}


@dataclass(frozen=True)
class Config:
    project_root: Path
    ff_dir: Path
    n_replica: int
    omp_threads: int
    native_restraint_strength: float
    rmsd_k: int
    minibatch_size: int
    n_frame: float
    sim_time: float
    alpha: float
    restart_from_last: bool
    worker_launch: str
    worker_python: str
    max_parallel_workers: int
    symlay_dense_grid_size: int
    symlay_support_margin: float
    symlay_projection_enabled: bool
    pair_ff_itp: Path
    pair_reg_enabled: bool
    pair_reg_cb_weight: float
    pair_reg_donor_weight: float
    pair_reg_acceptor_weight: float


@dataclass(frozen=True)
class MembraneNodeSpec:
    cb_node_name: str
    hb_node_name: str
    cb_coeff_shape: Tuple[int, ...]
    cb_inner_shape: Tuple[int, ...]
    hb_coeff_shape: Tuple[int, ...]
    hb_inner_shape: Tuple[int, ...]
    cb_coeff_size: int
    cb_inner_size: int
    hb_coeff_size: int
    hb_inner_size: int
    cb_active_size: int
    hb_active_size: int


def build_config() -> Config:
    ff_dir = Path(
        os.environ.get("CONDIV_FF_DIR", str(PROJECT_ROOT / "parameters" / "ff_2.1"))
    ).expanduser().resolve()
    worker_launch = os.environ.get("CONDIV_WORKER_LAUNCH", "auto").strip().lower()
    if worker_launch not in {"auto", "local", "srun"}:
        raise ValueError("CONDIV_WORKER_LAUNCH must be one of auto|local|srun")
    legacy_threads = int(os.environ.get("CONDIV_N_THREADS", "8"))
    slurm_active = bool(os.environ.get("SLURM_JOB_ID"))
    n_replica_default = "8" if slurm_active else str(legacy_threads)
    n_replica = int(os.environ.get("CONDIV_N_REPLICA", n_replica_default))
    if slurm_active:
        omp_threads_default = str(int(os.environ.get("SLURM_CPUS_PER_TASK", "1")))
    else:
        omp_threads_default = str(legacy_threads if "CONDIV_N_REPLICA" not in os.environ else 1)

    return Config(
        project_root=PROJECT_ROOT,
        ff_dir=ff_dir,
        n_replica=n_replica,
        omp_threads=int(os.environ.get("CONDIV_OMP_THREADS", omp_threads_default)),
        native_restraint_strength=float(os.environ.get("CONDIV_NATIVE_RESTRAINT_STRENGTH", str(1.0 / 3.0**2))),
        rmsd_k=int(os.environ.get("CONDIV_RMSD_K", "15")),
        minibatch_size=int(os.environ.get("CONDIV_MINIBATCH_SIZE", "15")),
        n_frame=float(os.environ.get("CONDIV_N_FRAME", "200")),
        sim_time=float(os.environ.get("CONDIV_SIM_TIME", "2000")),
        alpha=float(os.environ.get("CONDIV_ALPHA", "1.0")),
        restart_from_last=_env_bool("CONDIV_RESTART_FROM_LAST", False),
        worker_launch=worker_launch,
        worker_python=os.environ.get("CONDIV_WORKER_PYTHON", sys.executable),
        max_parallel_workers=int(os.environ.get("CONDIV_MAX_PARALLEL_WORKERS", "0")),
        symlay_dense_grid_size=int(
            os.environ.get("CONDIV_SYMLAY_DENSE_GRID_SIZE", str(DEFAULT_DENSE_GRID_SIZE))
        ),
        symlay_support_margin=float(
            os.environ.get("CONDIV_SYMLAY_SUPPORT_MARGIN", str(DEFAULT_SUPPORT_MARGIN))
        ),
        symlay_projection_enabled=_env_bool("CONDIV_SYMLAY_PROJECTION_ENABLED", True),
        pair_ff_itp=Path(
            os.environ.get(
                "CONDIV_PAIR_FF_ITP",
                str(PROJECT_ROOT / "example" / "16.MARTINI" / "ff_dry" / "dry_martini_v2.1.itp"),
            )
        ).expanduser().resolve(),
        pair_reg_enabled=_env_bool("CONDIV_PAIR_REG_ENABLED", True),
        pair_reg_cb_weight=float(os.environ.get("CONDIV_PAIR_REG_CB_WEIGHT", "0.05")),
        pair_reg_donor_weight=float(os.environ.get("CONDIV_PAIR_REG_DONOR_WEIGHT", "0.10")),
        pair_reg_acceptor_weight=float(os.environ.get("CONDIV_PAIR_REG_ACCEPTOR_WEIGHT", "0.10")),
    )


CONFIG = build_config()


resnames = [
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
]


Target = collections.namedtuple(
    "Target", "fasta native native_path breakfile_path thickness nail_file init_path n_res chi".split()
)
UpdateBase = collections.namedtuple("UpdateBase", "cb icb hb ihb".split())


class Update(UpdateBase):
    def __init__(self, *args):
        super().__init__()

    def _do_binary(self, other, op):
        try:
            len(other)
            is_seq = True
        except TypeError:
            is_seq = False

        if is_seq:
            assert len(self) == len(other)
            ret = []
            for a, b in zip(self, other):
                if a is None or b is None:
                    ret.append(None)
                else:
                    ret.append(op(a, b))
        else:
            ret = [None if a is None else op(a, other) for a in self]
        return Update(*ret)

    def __add__(self, other):
        return self._do_binary(other, lambda a, b: a + b)

    def __sub__(self, other):
        return self._do_binary(other, lambda a, b: a - b)

    def __mul__(self, other):
        return self._do_binary(other, lambda a, b: a * b)

    def __truediv__(self, other):
        return self._do_binary(other, lambda a, b: a / b)


class AdamSolver:
    """Simple Adam solver for membrane parameter blocks."""

    def __init__(self, n_comp: int, alpha=1e-2, beta1=0.8, beta2=0.96, epsilon=1e-6):
        self.n_comp = n_comp
        self.alpha = alpha
        self.beta1 = beta1
        self.beta2 = beta2
        self.epsilon = epsilon
        self.step_num = 0
        self.grad1 = [0.0 for _ in range(n_comp)]
        self.grad2 = [0.0 for _ in range(n_comp)]

    @staticmethod
    def _read_comp(value, idx):
        try:
            return value[idx]
        except Exception:
            return value

    def update_for_d_obj(self):
        return [0.0 for _ in self.grad1]

    def update_step(self, grad):
        self.step_num += 1
        t = self.step_num
        update = [None] * len(self.grad1)

        for i, g in enumerate(grad):
            b1 = self._read_comp(self.beta1, i)
            b2 = self._read_comp(self.beta2, i)
            a = self._read_comp(self.alpha, i)
            eps = self._read_comp(self.epsilon, i)

            self.grad1[i] = b1 * self.grad1[i] + (1.0 - b1) * g
            self.grad2[i] = b2 * self.grad2[i] + (1.0 - b2) * g**2

            grad1corr = self.grad1[i] / (1 - b1**t)
            grad2corr = self.grad2[i] / (1 - b2**t)
            update[i] = -a * grad1corr / (np.sqrt(grad2corr) + eps)

        return update

    def log_state(self, direc: str):
        with open(os.path.join(direc, "solver_state.pkl"), "wb") as fh:
            cp.dump(
                dict(
                    step_num=self.step_num,
                    grad1=self.grad1,
                    grad2=self.grad2,
                    solver=str(self),
                ),
                fh,
                -1,
            )

    def __repr__(self):
        return (
            f"AdamSolver({self.n_comp}, alpha={self.alpha!r}, beta1={self.beta1!r}, "
            f"beta2={self.beta2!r}, epsilon={self.epsilon!r})"
        )


# Stabilize pickling for checkpoint portability.
Target.__module__ = "ConDiv_mem"
Update.__module__ = "ConDiv_mem"
AdamSolver.__module__ = "ConDiv_mem"


def save_checkpoint(path: str, state: dict) -> None:
    with open(path, "wb") as fh:
        cp.dump(state, fh, -1)


def load_checkpoint(path: str) -> dict:
    with open(path, "rb") as fh:
        return cp.load(fh)


def _coerce_membrane_file(init_dir: str) -> str:
    p = Path(init_dir).expanduser().resolve()
    if p.is_file():
        return str(p)
    return str((p / "membrane.h5").resolve())


def get_init_param(init_dir: str):
    init_param_files = dict(memb=_coerce_membrane_file(init_dir))
    with tb.open_file(init_param_files["memb"]) as t:
        cb_memb = t.root.cb_energy[:]
        icb_memb = t.root.icb_energy[:]
        hb_memb = t.root.hb_energy[:]
        ihb_memb = t.root.ihb_energy[:]

    param = Update(*([None] * 4))._replace(cb=cb_memb, icb=icb_memb, hb=hb_memb, ihb=ihb_memb)
    return param, init_param_files


def _resolve_layer_manifest_path(base_dir: str) -> Path:
    raw = os.environ.get("CONDIV_SYMLAY_LAYER_MANIFEST", "")
    if raw.strip():
        return Path(raw).expanduser().resolve()
    return (Path(base_dir).expanduser().resolve() / "layer_manifest.json").resolve()


def _project_update_to_symlay(
    param: Update,
    init_membrane_file: str,
    layer_manifest: dict,
    dense_grid_size: int,
    support_margin: float,
) -> Tuple[Update, dict]:
    arrays = load_membrane_arrays(Path(init_membrane_file))
    arrays["cb"] = np.asarray(param.cb, dtype=np.float64)
    arrays["icb"] = np.asarray(param.icb, dtype=np.float64)
    arrays["hb"] = np.asarray(param.hb, dtype=np.float64)
    arrays["ihb"] = np.asarray(param.ihb, dtype=np.float64)
    projected, summary = project_membrane_arrays(
        arrays,
        layer_manifest,
        dense_grid_size=dense_grid_size,
        support_margin=support_margin,
    )
    projected_param = Update(
        np.asarray(projected["cb"], dtype=np.float32),
        np.asarray(projected["icb"], dtype=np.float32),
        np.asarray(projected["hb"], dtype=np.float32),
        np.asarray(projected["ihb"], dtype=np.float32),
    )
    projection_stats = {
        "cb_max_projection_residual": float(
            max((row["projection_residual_max_abs"] for row in summary["cb"]), default=0.0)
        ),
        "icb_max_projection_residual": float(
            max((row["projection_residual_max_abs"] for row in summary["icb"]), default=0.0)
        ),
        "hb_max_projection_residual": float(
            max((row["projection_residual_max_abs"] for row in summary["hb"]), default=0.0)
        ),
        "ihb_max_projection_residual": float(
            max((row["projection_residual_max_abs"] for row in summary["ihb"]), default=0.0)
        ),
        "cb_max_symmetry_residual": float(
            max((row["symmetry_residual_max_abs"] for row in summary["cb"]), default=0.0)
        ),
        "icb_max_symmetry_residual": float(
            max((row["symmetry_residual_max_abs"] for row in summary["icb"]), default=0.0)
        ),
        "hb_max_symmetry_residual": float(
            max((row["symmetry_residual_max_abs"] for row in summary["hb"]), default=0.0)
        ),
        "ihb_max_symmetry_residual": float(
            max((row["symmetry_residual_max_abs"] for row in summary["ihb"]), default=0.0)
        ),
    }
    return projected_param, projection_stats


def _seed_constrained_membrane(
    source_init: str,
    base_dir: str,
    layer_manifest: dict,
    dense_grid_size: int,
    support_margin: float,
) -> Tuple[str, dict]:
    source_membrane = _coerce_membrane_file(source_init)
    seed_dir = Path(base_dir).expanduser().resolve() / "seed_forcefield"
    seed_dir.mkdir(parents=True, exist_ok=True)
    seed_file = seed_dir / "membrane.h5"

    arrays = load_membrane_arrays(Path(source_membrane))
    projected, summary = project_membrane_arrays(
        arrays,
        layer_manifest,
        dense_grid_size=dense_grid_size,
        support_margin=support_margin,
    )
    write_projected_membrane_file(Path(source_membrane), seed_file, projected)

    seed_summary = {
        "schema": "condiv_symlay_seed_summary_v1",
        "source_membrane": repo_relative(Path(source_membrane)),
        "seed_membrane": repo_relative(seed_file),
        "layer_manifest": repo_relative(_resolve_layer_manifest_path(base_dir)),
        "cb_support": [float(projected["cb_z_min"]), float(projected["cb_z_max"])],
        "hb_support": [float(projected["hb_z_min"]), float(projected["hb_z_max"])],
        "cb_max_projection_residual": float(
            max((row["projection_residual_max_abs"] for row in summary["cb"]), default=0.0)
        ),
        "hb_max_projection_residual": float(
            max((row["projection_residual_max_abs"] for row in summary["hb"]), default=0.0)
        ),
        "cb_max_symmetry_residual": float(
            max((row["symmetry_residual_max_abs"] for row in summary["cb"]), default=0.0)
        ),
        "hb_max_symmetry_residual": float(
            max((row["symmetry_residual_max_abs"] for row in summary["hb"]), default=0.0)
        ),
        "dense_grid_size": int(dense_grid_size),
        "support_margin": float(support_margin),
    }
    (seed_dir / "seed_summary.json").write_text(
        json.dumps(seed_summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return str(seed_file), seed_summary


def _shape_tuple(shape_like) -> Tuple[int, ...]:
    return tuple(int(x) for x in shape_like)


def _membrane_shapes_from_config(config_base: str) -> MembraneNodeSpec:
    with tb.open_file(config_base) as t:
        pot = t.root.input.potential
        if hasattr(pot, "cb_surf_membrane_potential"):
            cb_node_name = "cb_surf_membrane_potential"
            hb_node_name = "hb_surf_membrane_potential"
            cb_coeff_shape = _shape_tuple(pot.cb_surf_membrane_potential.coeff.shape)
            cb_inner_shape = _shape_tuple(pot.cb_surf_membrane_potential.coeff_inner.shape)
            hb_coeff_shape = _shape_tuple(pot.hb_surf_membrane_potential.coeff.shape)
            hb_inner_shape = _shape_tuple(pot.hb_surf_membrane_potential.coeff_inner.shape)
        else:
            cb_node_name = "cb_membrane_potential"
            hb_node_name = "hb_membrane_potential"
            cb_coeff_shape = _shape_tuple(pot.cb_membrane_potential.coeff.shape)
            cb_inner_shape = (0,)
            hb_coeff_shape = _shape_tuple(pot.hb_membrane_potential.coeff.shape)
            hb_inner_shape = (0,)

    cb_coeff_size = int(np.prod(cb_coeff_shape))
    cb_inner_size = int(np.prod(cb_inner_shape))
    hb_coeff_size = int(np.prod(hb_coeff_shape))
    hb_inner_size = int(np.prod(hb_inner_shape))

    return MembraneNodeSpec(
        cb_node_name=cb_node_name,
        hb_node_name=hb_node_name,
        cb_coeff_shape=cb_coeff_shape,
        cb_inner_shape=cb_inner_shape,
        hb_coeff_shape=hb_coeff_shape,
        hb_inner_shape=hb_inner_shape,
        cb_coeff_size=cb_coeff_size,
        cb_inner_size=cb_inner_size,
        hb_coeff_size=hb_coeff_size,
        hb_inner_size=hb_inner_size,
        cb_active_size=cb_coeff_size,
        hb_active_size=hb_coeff_size,
    )


def _resolve_active_param_size(engine: ue.Upside, node_name: str, coeff_size: int, inner_size: int) -> int:
    sizes_to_try = [coeff_size]
    total_size = coeff_size + inner_size
    if inner_size and total_size not in sizes_to_try:
        sizes_to_try.append(total_size)

    errors = []
    for size in sizes_to_try:
        try:
            engine.get_param((size,), node_name)
            return int(size)
        except RuntimeError as exc:
            errors.append(f"{size}: {exc}")

    raise RuntimeError(
        f"Unable to resolve active parameter size for node {node_name}. "
        f"Tried {sizes_to_try}. Errors: {' | '.join(errors)}"
    )


def resolve_membrane_param_spec(config_base: str) -> MembraneNodeSpec:
    raw = _membrane_shapes_from_config(config_base)
    engine = ue.Upside(config_base)
    cb_active = _resolve_active_param_size(engine, raw.cb_node_name, raw.cb_coeff_size, raw.cb_inner_size)
    hb_active = _resolve_active_param_size(engine, raw.hb_node_name, raw.hb_coeff_size, raw.hb_inner_size)

    valid_cb_sizes = {raw.cb_coeff_size, raw.cb_coeff_size + raw.cb_inner_size}
    valid_hb_sizes = {raw.hb_coeff_size, raw.hb_coeff_size + raw.hb_inner_size}
    if cb_active not in valid_cb_sizes:
        raise RuntimeError(
            f"Unexpected active size for {raw.cb_node_name}: {cb_active}; "
            f"expected one of {sorted(valid_cb_sizes)}"
        )
    if hb_active not in valid_hb_sizes:
        raise RuntimeError(
            f"Unexpected active size for {raw.hb_node_name}: {hb_active}; "
            f"expected one of {sorted(valid_hb_sizes)}"
        )

    return MembraneNodeSpec(
        cb_node_name=raw.cb_node_name,
        hb_node_name=raw.hb_node_name,
        cb_coeff_shape=raw.cb_coeff_shape,
        cb_inner_shape=raw.cb_inner_shape,
        hb_coeff_shape=raw.hb_coeff_shape,
        hb_inner_shape=raw.hb_inner_shape,
        cb_coeff_size=raw.cb_coeff_size,
        cb_inner_size=raw.cb_inner_size,
        hb_coeff_size=raw.hb_coeff_size,
        hb_inner_size=raw.hb_inner_size,
        cb_active_size=cb_active,
        hb_active_size=hb_active,
    )


def _validate_membrane_param_shapes(param: Update, membrane_file: str) -> None:
    with tb.open_file(membrane_file) as t:
        expected = {
            "cb": _shape_tuple(t.root.cb_energy.shape),
            "icb": _shape_tuple(t.root.icb_energy.shape),
            "hb": _shape_tuple(t.root.hb_energy.shape),
            "ihb": _shape_tuple(t.root.ihb_energy.shape),
        }
    actual = {
        "cb": _shape_tuple(param.cb.shape),
        "icb": _shape_tuple(param.icb.shape),
        "hb": _shape_tuple(param.hb.shape),
        "ihb": _shape_tuple(param.ihb.shape),
    }
    mismatched = [k for k in expected if expected[k] != actual[k]]
    if mismatched:
        details = "; ".join(f"{k}: expected {expected[k]}, got {actual[k]}" for k in mismatched)
        raise RuntimeError(f"Membrane parameter shape mismatch against force field: {details}")


def _validate_membrane_file_vs_model(config_base: str, membrane_file: str) -> MembraneNodeSpec:
    spec = resolve_membrane_param_spec(config_base)
    with tb.open_file(membrane_file) as t:
        ff_shapes = {
            "cb_coeff": _shape_tuple(t.root.cb_energy.shape),
            "cb_inner": _shape_tuple(t.root.icb_energy.shape),
            "hb_coeff": _shape_tuple(t.root.hb_energy.shape),
            "hb_inner": _shape_tuple(t.root.ihb_energy.shape),
        }

    expected_shapes = {
        "cb_coeff": spec.cb_coeff_shape,
        "hb_coeff": spec.hb_coeff_shape,
    }
    mismatched = [k for k in expected_shapes if ff_shapes[k] != expected_shapes[k]]
    if spec.cb_inner_size > 0 and ff_shapes["cb_inner"] != spec.cb_inner_shape:
        mismatched.append("cb_inner")
    if spec.hb_inner_size > 0 and ff_shapes["hb_inner"] != spec.hb_inner_shape:
        mismatched.append("hb_inner")
    if mismatched:
        details = "; ".join(
            f"{k}: ff {ff_shapes[k]} vs model {expected_shapes.get(k, spec.cb_inner_shape if k == 'cb_inner' else spec.hb_inner_shape)}"
            for k in mismatched
        )
        raise RuntimeError(f"Membrane force-field/model dimension mismatch: {details}")

    return spec


def shift_inner_curve(param_inner: np.ndarray, param_outer: np.ndarray) -> np.ndarray:
    spoint = 8
    n_node = param_inner.size
    xx = np.arange(1, n_node - 1)
    y_i = ue.clamped_spline_value(param_inner, xx)
    y_o = ue.clamped_spline_value(param_outer, xx)
    y_i[spoint:] += y_o[spoint] - y_i[spoint]
    y_i[:spoint] = y_o[:spoint]
    return ue.clamped_spline_solve(y_i)


def expand_param(params: Update, orig_param_files: Dict[str, str], new_param_files: Dict[str, str]) -> None:
    shutil.copyfile(orig_param_files["memb"], new_param_files["memb"])

    n_aa = params.cb.shape[0]
    params_icb = np.zeros_like(params.icb)
    for aa in range(n_aa):
        params_icb[aa] = shift_inner_curve(params.icb[aa], params.cb[aa, 0])

    n_tp = params.hb.shape[0]
    n_op = params.hb.shape[1]
    params_ihb = np.zeros_like(params.ihb)
    for tp in range(n_tp):
        for op in range(n_op):
            params_ihb[tp, op] = shift_inner_curve(params.ihb[tp, op], params.hb[tp, op])

    manifest_path = os.environ.get("CONDIV_SYMLAY_LAYER_MANIFEST", "").strip()
    if CONFIG.symlay_projection_enabled and manifest_path:
        layer_manifest = load_layer_manifest(Path(manifest_path))
        arrays = load_membrane_arrays(Path(orig_param_files["memb"]))
        arrays["cb"] = np.asarray(params.cb, dtype=np.float64)
        arrays["icb"] = np.asarray(params_icb, dtype=np.float64)
        arrays["hb"] = np.asarray(params.hb, dtype=np.float64)
        arrays["ihb"] = np.asarray(params_ihb, dtype=np.float64)
        projected, _ = project_membrane_arrays(
            arrays,
            layer_manifest,
            dense_grid_size=CONFIG.symlay_dense_grid_size,
            support_margin=CONFIG.symlay_support_margin,
        )
        cb_energy = projected["cb"]
        icb_energy = projected["icb"]
        hb_energy = projected["hb"]
        ihb_energy = projected["ihb"]
    else:
        cb_energy = np.asarray(params.cb, dtype=np.float64)
        icb_energy = np.asarray(params_icb, dtype=np.float64)
        hb_energy = np.asarray(params.hb, dtype=np.float64)
        ihb_energy = np.asarray(params_ihb, dtype=np.float64)

    with tb.open_file(new_param_files["memb"], "a") as t:
        t.root.cb_energy[:] = np.asarray(cb_energy, dtype=np.float32)
        t.root.icb_energy[:] = np.asarray(icb_energy, dtype=np.float32)
        t.root.hb_energy[:] = np.asarray(hb_energy, dtype=np.float32)
        t.root.ihb_energy[:] = np.asarray(ihb_energy, dtype=np.float32)


def gen_swap_set(st: int, n_rep: int) -> Tuple[str, str]:
    a = ""
    for i in range(st, n_rep - 1, 2):
        a += f"{i}-{i + 1},"
    b = ""
    for i in range(st + 1, n_rep - 1, 2):
        b += f"{i}-{i + 1},"
    return a[:-1], b[:-1]


def _temperature_schedule(n_rep: int) -> np.ndarray:
    if n_rep <= 1:
        return np.array([0.8, 0.8], dtype=float)
    if n_rep == 8:
        return np.array([0.8, 0.8, 0.82, 0.84, 0.86, 0.88, 0.90, 0.92], dtype=float)
    arr = 0.8 + 0.12 * np.arange(n_rep, dtype=float) / max(1, n_rep - 1)
    arr[0] = 0.8
    if n_rep > 1:
        arr[1] = 0.8
    return arr


def _forcefield_kwargs(ff_dir: str, membrane_file: str, thickness: str, chain_break: str) -> dict:
    ff = Path(ff_dir)
    common = CONFIG.project_root / "parameters" / "common"
    kwargs = dict(
        environment_potential=str(ff / "environment.h5"),
        environment_potential_type=1,
        bb_environment_potential=str(ff / "bb_env.dat"),
        rotamer_interaction=str(ff / "sidechain.h5"),
        rotamer_placement=str(ff / "sidechain.h5"),
        dynamic_rotamer_1body=True,
        hbond_energy=str(ff / "hbond.h5"),
        rama_sheet_mix_energy=str(ff / "sheet"),
        rama_param_deriv=True,
        rama_library=str(common / "rama.dat"),
        reference_state_rama=str(common / "rama_reference.pkl"),
        membrane_potential=membrane_file,
        membrane_thickness=float(thickness),
        surface=True,
    )
    if chain_break and os.path.exists(chain_break):
        kwargs["chain_break_from_file"] = chain_break
    return kwargs


def _prepare_chain_break_file(chain_break_path: str, output_path: str) -> str:
    if not chain_break_path or not os.path.exists(chain_break_path):
        return ""

    with open(chain_break_path, "r", encoding="utf-8") as fh:
        lines = [line.strip() for line in fh if line.strip()]
    if not lines:
        return ""
    if len(lines) >= 2:
        return chain_break_path

    chain_first_residue = [int(x) for x in lines[0].split()]
    chain_counts = ["1"] * (len(chain_first_residue) + 1)
    with open(output_path, "w", encoding="utf-8") as fh:
        fh.write(" ".join(str(x) for x in chain_first_residue) + "\n")
        fh.write(" ".join(chain_counts) + "\n")
    return output_path


def _fallback_inner_shapes_from_config(config_base: str) -> Tuple[Tuple[int, ...], Tuple[int, ...]]:
    default_shapes = ((0,), (0,))
    try:
        with tb.open_file(config_base) as t:
            args = t.root.input.args._v_attrs
            membrane_file = getattr(args, "channel_membrane_potential", "") or getattr(args, "membrane_potential", "")
        if not membrane_file or not os.path.exists(membrane_file):
            return default_shapes
        with tb.open_file(membrane_file) as mt:
            return _shape_tuple(mt.root.icb_energy.shape), _shape_tuple(mt.root.ihb_energy.shape)
    except Exception:
        return default_shapes


def compute_divergence(config_base: str, pos: np.ndarray, mode: int = 0) -> Update:
    del mode
    spec = resolve_membrane_param_spec(config_base)
    cb_inner_shape = spec.cb_inner_shape
    hb_inner_shape = spec.hb_inner_shape
    if spec.cb_inner_size == 0 or spec.hb_inner_size == 0:
        fb_cb_inner_shape, fb_hb_inner_shape = _fallback_inner_shapes_from_config(config_base)
        if spec.cb_inner_size == 0:
            cb_inner_shape = fb_cb_inner_shape
        if spec.hb_inner_size == 0:
            hb_inner_shape = fb_hb_inner_shape

    engine = ue.Upside(config_base)
    contrast = Update([], [], [], [])

    for i in range(pos.shape[0]):
        engine.energy(pos[i])

        dp_cb_memb = engine.get_param_deriv((spec.cb_active_size,), spec.cb_node_name)
        if spec.cb_active_size == spec.cb_coeff_size:
            odp_cb_memb = dp_cb_memb.reshape(spec.cb_coeff_shape)
            idp_cb_memb = np.zeros(cb_inner_shape, dtype=np.float32)
        else:
            odp_cb_memb = dp_cb_memb[: spec.cb_coeff_size].reshape(spec.cb_coeff_shape)
            idp_cb_memb = dp_cb_memb[spec.cb_coeff_size :].reshape(spec.cb_inner_shape)

        dp_hb_memb = engine.get_param_deriv((spec.hb_active_size,), spec.hb_node_name)
        if spec.hb_active_size == spec.hb_coeff_size:
            odp_hb_memb = dp_hb_memb.reshape(spec.hb_coeff_shape)
            idp_hb_memb = np.zeros(hb_inner_shape, dtype=np.float32)
        else:
            odp_hb_memb = dp_hb_memb[: spec.hb_coeff_size].reshape(spec.hb_coeff_shape)
            idp_hb_memb = dp_hb_memb[spec.hb_coeff_size :].reshape(spec.hb_inner_shape)
        if not np.isfinite(odp_cb_memb).all() or not np.isfinite(idp_cb_memb).all():
            raise RuntimeError(f"NONFINITE_MEMBRANE_DERIV {spec.cb_node_name}")
        if not np.isfinite(odp_hb_memb).all() or not np.isfinite(idp_hb_memb).all():
            raise RuntimeError(f"NONFINITE_MEMBRANE_DERIV {spec.hb_node_name}")

        contrast.cb.append(odp_cb_memb)
        contrast.icb.append(idp_cb_memb)
        contrast.hb.append(odp_hb_memb)
        contrast.ihb.append(idp_hb_memb)

    return Update(*[np.asarray(x) for x in contrast])


def _safe_energy_component(engine: ue.Upside, node_name: str) -> float:
    try:
        return float(engine.get_output(node_name)[0, 0])
    except Exception:
        return float("nan")


def _select_finite_equil_frames(pos: np.ndarray, equil_fraction: float) -> np.ndarray:
    if pos.ndim != 3:
        raise RuntimeError(f"Unexpected trajectory shape {pos.shape}, expected (n_frame,n_atom,3)")
    n_frame = int(pos.shape[0])
    if n_frame <= 0:
        raise RuntimeError("No output frames available")

    start = min(int(equil_fraction * n_frame), max(0, n_frame - 1))
    sliced = pos[start:]
    mask = np.isfinite(sliced).all(axis=(1, 2))
    filtered = sliced[mask]
    if filtered.shape[0] > 0:
        return filtered

    mask_all = np.isfinite(pos).all(axis=(1, 2))
    filtered_all = pos[mask_all]
    if filtered_all.shape[0] > 0:
        return filtered_all

    raise RuntimeError("No finite trajectory frames available")


def _rmsd_atom_slice(n_atom: int, k: int) -> slice:
    k_max = max(0, n_atom // 2 - 1)
    k_eff = max(0, min(k, k_max))
    if k_eff == 0:
        return slice(None)
    return slice(k_eff, -k_eff)


def compute_frame_properties(argv):
    path_code, config_base, trajs, rep_id, start, h_thickness, rmsd_k = argv
    del path_code
    del rep_id

    fn_ref = mu.load_upside_ref(config_base, add_atoms=False)
    top = fn_ref.top
    ca = top.select("name CA")
    fn_ref_ca = fn_ref.atom_slice(ca)

    z_ref = fn_ref.xyz[0][:, 2] * 10.0
    ndx_tm = np.where((z_ref < h_thickness) * (z_ref > -h_thickness))[0]
    sz_ref = z_ref[ndx_tm]
    a1, b1 = np.histogram(sz_ref, 20, density=True)
    c1 = (b1[1:] + b1[:-1]) * 0.5

    fn = mu.load_upside_traj(trajs, add_atoms=False)
    xyz = fn.xyz[:] * 10.0
    fn_ca = fn.atom_slice(ca)
    sn = fn.n_frames
    rmsd = ru.traj_rmsd(fn_ca.xyz[:, rmsd_k:-rmsd_k, :] * 10.0, fn_ref_ca.xyz[0, rmsd_k:-rmsd_k, :] * 10.0)
    rg = md.compute_rg(fn_ca)
    del rg

    z_traj = xyz[start:, :, 2]
    sz_traj = z_traj[:, ndx_tm]
    a2, b2 = np.histogram(sz_traj.flatten(), 20, density=True)
    c2 = (b2[1:] + b2[:-1]) * 0.5

    pot = [
        "rotamer",
        "sigmoid_coupling_environment",
        "bb_sigmoid_coupling_environment",
        "hbond_energy",
        "rama_map_pot",
        "cb_surf_membrane_potential",
        "hb_surf_membrane_potential",
    ]
    energies = np.zeros((len(pot) + 1, sn), dtype=np.float64)
    engine = ue.Upside(config_base)
    for i in range(sn):
        energies[0, i] = engine.energy(xyz[i])
        for j, node_name in enumerate(pot):
            energies[j + 1, i] = _safe_energy_component(engine, node_name)

    return [rmsd, energies, c1, a1, c2, a2]


def compute_frame_divergence(argv):
    config_base, trajs, mode, start = argv
    with tb.open_file(trajs) as t:
        xyz = t.root.output.pos[start:, 0]
    div = compute_divergence(config_base, xyz, mode)
    return [x for x in div]


def _has_output_frames(config_path: str) -> bool:
    try:
        with tb.open_file(config_path) as t:
            if not hasattr(t.root, "output"):
                return False
            return int(t.root.output.pos.shape[0]) > 0
    except Exception:
        return False


def _choose_worker_launch(mode: str) -> str:
    if mode == "auto":
        return "srun" if os.environ.get("SLURM_JOB_ID") else "local"
    return mode


def _apply_runtime_config(state: dict) -> dict:
    state["project_root"] = str(CONFIG.project_root)
    state["ff_dir"] = str(CONFIG.ff_dir)
    state["worker_launch"] = CONFIG.worker_launch
    state["worker_python"] = CONFIG.worker_python
    state["n_replica"] = CONFIG.n_replica
    state["omp_threads"] = CONFIG.omp_threads
    state["max_parallel_workers"] = CONFIG.max_parallel_workers
    if __file__:
        state["worker_path"] = str(Path(__file__).resolve())
    state["symlay_dense_grid_size"] = CONFIG.symlay_dense_grid_size
    state["symlay_support_margin"] = CONFIG.symlay_support_margin
    state["symlay_projection_enabled"] = CONFIG.symlay_projection_enabled
    state["pair_ff_itp"] = str(CONFIG.pair_ff_itp)
    state["pair_reg_enabled"] = CONFIG.pair_reg_enabled
    state["pair_reg_cb_weight"] = CONFIG.pair_reg_cb_weight
    state["pair_reg_donor_weight"] = CONFIG.pair_reg_donor_weight
    state["pair_reg_acceptor_weight"] = CONFIG.pair_reg_acceptor_weight
    if state.get("pair_reg_enabled", False) and "pair_teacher" not in state:
        state["pair_teacher"] = build_pair_teacher(state["layer_manifest"], Path(state["pair_ff_itp"]))
    return state


def _run_worker_subprocess(
    state: dict,
    nm: str,
    direc: str,
    worker_argv: List[str],
) -> Tuple[sp.Popen, Optional[object]]:
    # Keep the Python worker local. When Slurm is active, the worker launches
    # the actual Upside replica bundle through `srun`, which lets one worker
    # reserve `n_replica` task slots like the replica-exchange reference flow.
    outfile = open(f"{direc}/{nm}.output_worker", "w", encoding="utf-8")
    cmd = [state["worker_python"], *worker_argv]
    return sp.Popen(cmd, close_fds=True, stdout=outfile, stderr=sp.STDOUT), outfile


def run_minibatch(
    state: dict,
    param: Update,
    initial_param_files: Dict[str, str],
    direc: str,
    minibatch,
    solver: AdamSolver,
    sim_time: float,
):
    os.makedirs(direc, exist_ok=True)
    print(direc)
    print()

    d_obj_param_files = {
        k: os.path.join(direc, f"nesterov_temp__{os.path.basename(x)}") for k, x in initial_param_files.items()
    }
    expand_param(param + solver.update_for_d_obj(), initial_param_files, d_obj_param_files)

    with open(os.path.join(direc, "sim_time"), "w", encoding="utf-8") as fh:
        print(sim_time, file=fh)

    rmsd = {}
    com = {}
    change = []

    max_parallel_workers = state.get("max_parallel_workers", 0)
    if max_parallel_workers <= 0:
        max_parallel_workers = len(minibatch)
    max_parallel_workers = min(max_parallel_workers, len(minibatch))

    print(
        "Worker dispatch=%s, replicas/worker=%i, omp_threads/upside=%i, max_parallel_workers=%i, minibatch_targets=%i"
        % (
            _choose_worker_launch(state["worker_launch"]),
            int(state["n_replica"]),
            int(state["omp_threads"]),
            max_parallel_workers,
            len(minibatch),
        )
    )

    for chunk_start in range(0, len(minibatch), max_parallel_workers):
        chunk = minibatch[chunk_start : chunk_start + max_parallel_workers]
        jobs = collections.OrderedDict()

        for nm, target in chunk[::-1]:
            worker_argv = [
                state["worker_path"],
                "worker",
                nm,
                direc,
                target.fasta,
                target.native_path,
                target.breakfile_path,
                str(target.thickness),
                target.nail_file,
                target.init_path,
                str(target.n_res),
                target.chi,
                cp.dumps(d_obj_param_files, protocol=4).hex(),
                str(sim_time),
                str(state["n_replica"]),
                str(state["omp_threads"]),
                _choose_worker_launch(state["worker_launch"]),
                state["project_root"],
                state["ff_dir"],
            ]
            jobs[nm] = _run_worker_subprocess(state, nm, direc, worker_argv)

        for nm, (job, outfile) in jobs.items():
            try:
                if job.wait() != 0:
                    print(nm, "WORKER_FAIL")
                    continue

                with open(f"{direc}/{nm}.divergence.pkl", "rb") as fh:
                    divergence = cp.load(fh)

                metrics = np.array(
                    [
                        float(divergence["rmsd_restrain"]),
                        float(divergence["rmsd"]),
                        float(divergence["com_restrain"]),
                        float(divergence["com"]),
                    ],
                    dtype=np.float64,
                )
                contrast = divergence["contrast"]
                contrast_finite = all(np.isfinite(np.asarray(x)).all() for x in contrast)
                if (not np.isfinite(metrics).all()) or (not contrast_finite):
                    print(nm, "NONFINITE_DIVERGENCE")
                    continue

                rmsd[nm] = (divergence["rmsd_restrain"], divergence["rmsd"])
                com[nm] = (divergence["com_restrain"], divergence["com"])
                change.append(contrast)
            finally:
                if outfile is not None:
                    outfile.close()

    if not change:
        raise RuntimeError("All jobs failed")

    with open(f"{direc}/rmsd.pkl", "wb") as fh:
        cp.dump(rmsd, fh, -1)

    rmsd_values = np.array(list(rmsd.values()))
    com_values = np.array(list(com.values()))
    print()
    print("Median RMSD %.2f %.2f" % tuple(np.median(rmsd_values, axis=0)))
    print("Median COM %.2f %.2f" % tuple(np.median(com_values, axis=0)))

    d_param = Update(*[np.sum(x, axis=0) for x in zip(*change)])
    pair_stats = {
        "pair_loss_cb": 0.0,
        "pair_loss_hb_donor": 0.0,
        "pair_loss_hb_acceptor": 0.0,
        "pair_loss_total": 0.0,
    }
    if state.get("pair_reg_enabled", False):
        d_pair, pair_stats = compute_pair_regularization_gradient(
            param,
            initial_param_files["memb"],
            state["layer_manifest"],
            state["pair_teacher"],
            float(state["pair_reg_cb_weight"]),
            float(state["pair_reg_donor_weight"]),
            float(state["pair_reg_acceptor_weight"]),
        )
        d_param = d_param + d_pair
    step_update = solver.update_step(d_param)
    grad_stats = {
        "n_success": int(len(change)),
        "n_total": int(len(minibatch)),
        "grad_norm_cb": float(np.linalg.norm(d_param.cb)),
        "grad_norm_icb": float(np.linalg.norm(d_param.icb)),
        "grad_norm_hb": float(np.linalg.norm(d_param.hb)),
        "grad_norm_ihb": float(np.linalg.norm(d_param.ihb)),
        "update_norm_cb": float(np.linalg.norm(step_update[0])),
        "update_norm_icb": float(np.linalg.norm(step_update[1])),
        "update_norm_hb": float(np.linalg.norm(step_update[2])),
        "update_norm_ihb": float(np.linalg.norm(step_update[3])),
        "pair_reg_enabled": bool(state.get("pair_reg_enabled", False)),
        "pair_teacher_schema": state.get("pair_teacher", {}).get("schema", ""),
        "pair_ff_itp": repo_relative(Path(state["pair_ff_itp"])) if state.get("pair_ff_itp") else "",
    }
    grad_stats.update(pair_stats)
    grad_stats["grad_norm_total"] = float(
        np.sqrt(
            grad_stats["grad_norm_cb"] ** 2
            + grad_stats["grad_norm_icb"] ** 2
            + grad_stats["grad_norm_hb"] ** 2
            + grad_stats["grad_norm_ihb"] ** 2
        )
    )
    grad_stats["update_norm_total"] = float(
        np.sqrt(
            grad_stats["update_norm_cb"] ** 2
            + grad_stats["update_norm_icb"] ** 2
            + grad_stats["update_norm_hb"] ** 2
            + grad_stats["update_norm_ihb"] ** 2
        )
    )

    new_param_files = {k: os.path.join(direc, os.path.basename(x)) for k, x in initial_param_files.items()}
    new_param = param + step_update
    projection_stats = {}
    if state.get("symlay_projection_enabled", True):
        new_param, projection_stats = _project_update_to_symlay(
            new_param,
            initial_param_files["memb"],
            state["layer_manifest"],
            int(state["symlay_dense_grid_size"]),
            float(state["symlay_support_margin"]),
        )
        grad_stats.update(projection_stats)
    expand_param(new_param, initial_param_files, new_param_files)
    solver.log_state(direc)

    with open(os.path.join(direc, "gradient_stats.json"), "w", encoding="utf-8") as fh:
        json.dump(grad_stats, fh, indent=2, sort_keys=True)

    return new_param, grad_stats


def main_worker():
    tstart = time.time()
    code = sys.argv[2]
    direc = sys.argv[3]
    fasta = sys.argv[4]
    native_path = sys.argv[5]
    chain_break = sys.argv[6]
    thickness = sys.argv[7]
    nail_file = sys.argv[8]
    init_path = sys.argv[9]
    n_res = int(sys.argv[10])
    chi = sys.argv[11]
    del nail_file
    del chi
    param_files = cp.loads(bytes.fromhex(sys.argv[12]))
    sim_time = float(sys.argv[13])
    n_replica = int(sys.argv[14])
    omp_threads = int(sys.argv[15])
    launch_mode = sys.argv[16]
    project_root = sys.argv[17]
    ff_dir = sys.argv[18]
    del project_root

    frame_interval = max(1, int(sim_time / CONFIG.n_frame))
    path_code = f"{direc}/{code}"

    chain_break_compat = _prepare_chain_break_file(
        chain_break, os.path.join(direc, f"{code}.chain_break_compat.txt")
    )
    kwargs = _forcefield_kwargs(ff_dir, param_files["memb"], thickness, chain_break_compat)
    T = _temperature_schedule(n_replica)

    try:
        config_base = f"{direc}/{code}.base.h5"
        configs = [re.sub(r"\.base\.h5", f".run.{i}.h5", config_base) for i in range(len(T))]

        initial_structure_npy = _materialize_initial_structure(
            native_path, os.path.join(direc, f"{code}.initial.npy")
        )
        ru.upside_config(fasta, config_base, initial_structure=initial_structure_npy, **kwargs)
        for i in range(1, len(T)):
            shutil.copyfile(config_base, configs[i])

        ru.upside_config(fasta, configs[0], initial_structure=initial_structure_npy, **kwargs)
        ru.advanced_config(
            configs[0],
            restraint_groups=[f"0-{n_res - 1}"],
            restraint_spring_constant=CONFIG.native_restraint_strength,
        )
        dim_spec = _validate_membrane_file_vs_model(configs[1], param_files["memb"])
        print(
            "Membrane dimensions "
            f"cb(coeff={dim_spec.cb_coeff_size}, inner={dim_spec.cb_inner_size}, active={dim_spec.cb_active_size}) "
            f"hb(coeff={dim_spec.hb_coeff_size}, inner={dim_spec.hb_inner_size}, active={dim_spec.hb_active_size})"
        )
    except Exception as exc:
        raise RuntimeError(f"CONFIG_FAIL: {exc}") from exc

    set1, set2 = gen_swap_set(1, len(T))
    swap_sets = [x for x in [set1, set2] if x]
    if not swap_sets:
        swap_sets = ru.swap_table2d(1, len(T))
    job = ru.run_upside(
        "srun" if launch_mode == "srun" else "",
        configs,
        sim_time,
        frame_interval,
        n_threads=omp_threads,
        temperature=T,
        swap_sets=swap_sets,
        mc_interval=5.0,
        replica_interval=5.0,
        time_step=0.01,
        disable_z_recentering=True,
        ntasks=n_replica,
    )
    run_retcode = job.job.wait()
    if run_retcode != 0:
        # Upside occasionally exits non-zero after writing usable output.
        if not all(_has_output_frames(cfg) for cfg in configs[: min(2, len(configs))]):
            raise RuntimeError("RUN_FAIL")
        print(f"WARNING: run_upside returned {run_retcode}, continuing with existing output frames")

    divergence = {}
    equil_fraction = 0.25

    with tb.open_file(configs[0]) as t:
        pos_restrain = _select_finite_equil_frames(t.root.output.pos[:, 0], equil_fraction)
    with tb.open_file(configs[1]) as t:
        pos_free = _select_finite_equil_frames(t.root.output.pos[:, 0], equil_fraction)

    target = _load_native_positions(initial_structure_npy)
    atom_slice = _rmsd_atom_slice(int(target.shape[0]), CONFIG.rmsd_k)
    divergence["rmsd_restrain"] = float(
        ru.traj_rmsd(pos_restrain[:, atom_slice], target[atom_slice]).mean()
    )
    divergence["rmsd"] = float(
        ru.traj_rmsd(pos_free[:, atom_slice], target[atom_slice]).mean()
    )
    divergence["com_restrain"] = float(pos_restrain[:, :, 2].mean())
    divergence["com"] = float(pos_free[:, :, 2].mean())

    all_pos = np.concatenate([pos_restrain, pos_free], axis=0)
    all_div = compute_divergence(configs[1], all_pos)
    n_pos0 = len(pos_restrain)
    divergence["contrast"] = Update(
        *[x[:n_pos0].mean(axis=0) - x[n_pos0:].mean(axis=0) for x in all_div]
    )
    metric_vals = np.array(
        [divergence["rmsd_restrain"], divergence["rmsd"], divergence["com_restrain"], divergence["com"]],
        dtype=np.float64,
    )
    if not np.isfinite(metric_vals).all():
        raise RuntimeError("NONFINITE_DIVERGENCE_METRICS")
    for term_name, term in zip(("cb", "icb", "hb", "ihb"), divergence["contrast"]):
        if not np.isfinite(term).all():
            raise RuntimeError(f"NONFINITE_DIVERGENCE_CONTRAST {term_name}")

    divergence["walltime"] = float(time.time() - tstart)

    if CONFIG.restart_from_last:
        with open(init_path, "rb") as ixyz:
            init_state = cp.load(ixyz)
        with open(init_path, "wb") as oxyz:
            for i in range(1, len(T)):
                with tb.open_file(configs[i]) as t:
                    init_state[i] = t.root.output.pos[-1:, 0]
            cp.dump(init_state, oxyz, -1)

    for fn in [config_base] + configs + [
        os.path.join(direc, f"{code}.initial.npy"),
        os.path.join(direc, f"{code}.chain_break_compat.txt"),
    ]:
        if os.path.exists(fn):
            os.remove(fn)

    with open(f"{direc}/{code}.divergence.pkl", "wb") as fh:
        cp.dump(divergence, fh, -1)


def _shuffle_minibatches(minibatches):
    all_items = list(itertools.chain.from_iterable(minibatches))
    np.random.shuffle(all_items)
    n_minibatch = len(minibatches)
    return [all_items[i::n_minibatch] for i in range(n_minibatch)]


def main_loop_iteration(state: dict) -> dict:
    print("#########################################")
    print("####      EPOCH %2i MINIBATCH %2i      ####" % (state["epoch"], state["i_mb"]))
    print("#########################################")
    print()
    sys.stdout.flush()

    tstart = time.time()
    state["mb_direc"] = os.path.join(state["base_dir"], f"epoch_{state['epoch']:02d}_minibatch_{state['i_mb']:02d}")
    if os.path.exists(state["mb_direc"]):
        shutil.rmtree(state["mb_direc"], ignore_errors=True)

    new_param, grad_stats = run_minibatch(
        state,
        state["param"],
        state["init_param_files"],
        state["mb_direc"],
        state["minibatches"][state["i_mb"]],
        state["solver"],
        state["sim_time"],
    )
    state["param"] = new_param
    state["last_gradient_stats"] = grad_stats

    print()
    print("%.0f seconds elapsed this minibatch" % (time.time() - tstart))
    print()
    sys.stdout.flush()

    state["i_mb"] += 1
    if state["i_mb"] >= len(state["minibatches"]):
        state["i_mb"] = 0
        state["epoch"] += 1
        state["minibatches"] = _shuffle_minibatches(state["minibatches"])

    return state


def main_loop(checkpoint_path: str, max_iter: int):
    checkpoint_path = os.path.abspath(checkpoint_path)
    for _ in range(max_iter):
        state = load_checkpoint(checkpoint_path)
        state = _apply_runtime_config(state)
        state = main_loop_iteration(state)
        checkpoint_path = os.path.join(state["mb_direc"], "checkpoint.pkl")
        save_checkpoint(checkpoint_path, state)


def _load_native_positions(path: str) -> np.ndarray:
    p = Path(path)
    if p.suffix == ".pkl":
        with open(path, "rb") as fh:
            arr = cp.load(fh, encoding="latin1")
    else:
        arr = np.load(path, allow_pickle=True)

    arr = np.asarray(arr)
    if arr.ndim == 3:
        arr = arr[:, :, 0]
    return arr.astype(np.float32, copy=False)


def _materialize_initial_structure(source_path: str, out_path: str) -> str:
    np.save(out_path, _load_native_positions(source_path))
    return out_path


def main_initialize(args):
    if len(args) != 4:
        raise RuntimeError("initialize requires: <init_dir> <protein_dir> <protein_list|cached> <base_dir>")

    state = {}
    state["init_dir"], state["protein_dir"], protein_list, state["base_dir"] = args
    os.makedirs(state["base_dir"], exist_ok=True)
    state["layer_manifest_path"] = str(_resolve_layer_manifest_path(state["base_dir"]))
    state["layer_manifest"] = load_layer_manifest(Path(state["layer_manifest_path"]))
    state["pair_ff_itp"] = str(CONFIG.pair_ff_itp)
    state["pair_reg_enabled"] = CONFIG.pair_reg_enabled
    state["pair_reg_cb_weight"] = CONFIG.pair_reg_cb_weight
    state["pair_reg_donor_weight"] = CONFIG.pair_reg_donor_weight
    state["pair_reg_acceptor_weight"] = CONFIG.pair_reg_acceptor_weight
    state["pair_teacher"] = build_pair_teacher(state["layer_manifest"], Path(state["pair_ff_itp"]))

    if not __file__:
        raise RuntimeError("No __file__ available")
    state["worker_path"] = os.path.join(state["base_dir"], "ConDiv_mem.py")
    state["worker_python"] = CONFIG.worker_python
    shutil.copy(__file__, state["worker_path"])
    shutil.copy(Path(__file__).with_name("symlay_utils.py"), os.path.join(state["base_dir"], "symlay_utils.py"))

    if protein_list != "cached":
        print("Reading training set")
        with open(protein_list, "r", encoding="utf-8") as fh:
            protein_names = [x.split()[0] for x in fh if x.strip()]
        assert protein_names[0] == "prot"
        protein_names = protein_names[1:]

        training_set = {}
        excluded = []
        for code in sorted(protein_names):
            base = os.path.join(state["protein_dir"], code)
            native_pos = _load_native_positions(base + ".initial.pkl")
            n_res = len(native_pos) // 3

            with open(base + ".thickness", "r", encoding="utf-8") as fh:
                thickness = fh.read().strip()

            max_sep = np.sqrt(np.sum(np.diff(native_pos, axis=0) ** 2, axis=-1)).max()
            if max_sep < 20000.0:
                break_path = base + ".chain_breaks"
                if not os.path.exists(break_path):
                    break_path = ""
                nail_path = base + ".nail"
                if not os.path.exists(nail_path):
                    nail_path = ""
                training_set[code] = Target(
                    base + ".fasta",
                    native_pos,
                    base + ".initial.pkl",
                    break_path,
                    thickness,
                    nail_path,
                    base + ".states.pkl",
                    n_res,
                    base + ".chi",
                )
            else:
                excluded.append(code)

        print("Excluded %i proteins due to chain breaks" % len(excluded))
        with open(os.path.join(state["base_dir"], "cd_training.pkl"), "wb") as fh:
            cp.dump(training_set, fh, -1)
    else:
        with open(os.path.join(state["base_dir"], "cd_training.pkl"), "rb") as fh:
            training_set = cp.load(fh)

    training_list = sorted(training_set.items(), key=lambda x: (x[1].n_res, x[0]))
    np.random.shuffle(training_list)

    minibatch_excess = len(training_list) % CONFIG.minibatch_size
    if minibatch_excess:
        training_list = training_list[:-minibatch_excess]
    n_minibatch = len(training_list) // CONFIG.minibatch_size
    if n_minibatch <= 0:
        raise RuntimeError("Not enough proteins to form one minibatch")

    state["minibatches"] = [
        training_list[i * CONFIG.minibatch_size : (i + 1) * CONFIG.minibatch_size] for i in range(n_minibatch)
    ]
    state["n_prot"] = n_minibatch * CONFIG.minibatch_size
    print(
        "Constructed %i minibatches of size %i (%i proteins)"
        % (n_minibatch, CONFIG.minibatch_size, state["n_prot"])
    )

    if state["init_dir"] != "cached":
        state["seed_source_init_dir"] = state["init_dir"]
        state["init_dir"], state["symlay_seed_summary"] = _seed_constrained_membrane(
            state["seed_source_init_dir"],
            state["base_dir"],
            state["layer_manifest"],
            CONFIG.symlay_dense_grid_size,
            CONFIG.symlay_support_margin,
        )
        state["param"], state["init_param_files"] = get_init_param(state["init_dir"])
        _validate_membrane_param_shapes(state["param"], state["init_param_files"]["memb"])
        with open(os.path.join(state["base_dir"], "condiv_init.pkl"), "wb") as fh:
            cp.dump(
                {
                    "init_dir": state["init_dir"],
                    "seed_source_init_dir": state["seed_source_init_dir"],
                    "param": state["param"],
                    "init_param_files": state["init_param_files"],
                    "layer_manifest_path": state["layer_manifest_path"],
                    "symlay_seed_summary": state["symlay_seed_summary"],
                    "pair_teacher": state["pair_teacher"],
                    "pair_ff_itp": state["pair_ff_itp"],
                },
                fh,
                -1,
            )
    else:
        with open(os.path.join(state["base_dir"], "condiv_init.pkl"), "rb") as fh:
            payload = cp.load(fh)
        if isinstance(payload, dict):
            state["init_dir"] = payload["init_dir"]
            state["param"] = payload["param"]
            state["init_param_files"] = payload["init_param_files"]
            state["seed_source_init_dir"] = payload.get("seed_source_init_dir", state["init_dir"])
            state["layer_manifest_path"] = payload.get("layer_manifest_path", state["layer_manifest_path"])
            state["symlay_seed_summary"] = payload.get("symlay_seed_summary", {})
            state["pair_teacher"] = payload.get("pair_teacher", state["pair_teacher"])
            state["pair_ff_itp"] = payload.get("pair_ff_itp", state["pair_ff_itp"])
        else:
            state["init_dir"], state["param"], state["init_param_files"] = payload
            state["seed_source_init_dir"] = state["init_dir"]
            state["symlay_seed_summary"] = {}
        state["layer_manifest"] = load_layer_manifest(Path(state["layer_manifest_path"]))
        if "pair_teacher" not in state or state["pair_teacher"].get("schema", "") != PAIR_TEACHER_SCHEMA:
            state["pair_teacher"] = build_pair_teacher(state["layer_manifest"], Path(state["pair_ff_itp"]))
        _validate_membrane_param_shapes(state["param"], state["init_param_files"]["memb"])

    state["initial_alpha"] = Update(0.02, 0.02, 0.02, 0.02) * CONFIG.alpha
    state["solver"] = AdamSolver(len(state["initial_alpha"]), alpha=state["initial_alpha"])
    state["sim_time"] = CONFIG.sim_time

    state = _apply_runtime_config(state)

    print()
    print("Optimizing with solver", state["solver"])
    print("ConDiv_symlay layer manifest", repo_relative(Path(state["layer_manifest_path"])))
    print("ConDiv_symlay slot sequence", "-".join(state["layer_manifest"]["full_type_sequence"]))
    print("ConDiv_symlay seed source", repo_relative(Path(state["seed_source_init_dir"])))
    print("ConDiv_symlay seed membrane", repo_relative(Path(state["init_param_files"]["memb"])))
    print("ConDiv_symlay pair teacher", state["pair_teacher"]["schema"], repo_relative(Path(state["pair_ff_itp"])))
    print(
        "ConDiv_symlay replica layout",
        f"replicas/worker={state['n_replica']}",
        f"omp_threads/upside={state['omp_threads']}",
    )
    print()

    state["epoch"] = 0
    state["i_mb"] = 0
    return state


def main_usage() -> str:
    return (
        "Usage:\n"
        "  ConDiv_mem.py initialize <init_dir> <protein_dir> <protein_list|cached> <base_dir>\n"
        "  ConDiv_mem.py restart <checkpoint.pkl> <max_iter>\n"
        "  ConDiv_mem.py worker <internal args...>\n"
    )


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise RuntimeError(main_usage())

    mode = sys.argv[1]
    if mode == "worker":
        main_worker()
    elif mode == "restart":
        if len(sys.argv[1:]) != 3:
            raise RuntimeError(main_usage())
        print("Running as PID %i on host %s" % (os.getpid(), socket.gethostname()))
        main_loop(sys.argv[2], int(sys.argv[3]))
    elif mode == "initialize":
        initial_state = main_initialize(sys.argv[2:])
        save_checkpoint(os.path.join(initial_state["base_dir"], "initial_checkpoint.pkl"), initial_state)
    else:
        raise RuntimeError(main_usage())
