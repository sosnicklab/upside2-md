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
    n_threads: int
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


def build_config() -> Config:
    ff_dir = Path(
        os.environ.get("CONDIV_FF_DIR", str(PROJECT_ROOT / "parameters" / "ff_2.1"))
    ).expanduser().resolve()
    worker_launch = os.environ.get("CONDIV_WORKER_LAUNCH", "auto").strip().lower()
    if worker_launch not in {"auto", "local", "srun"}:
        raise ValueError("CONDIV_WORKER_LAUNCH must be one of auto|local|srun")

    return Config(
        project_root=PROJECT_ROOT,
        ff_dir=ff_dir,
        n_threads=int(os.environ.get("CONDIV_N_THREADS", "8")),
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

    with tb.open_file(new_param_files["memb"], "a") as t:
        t.root.cb_energy[:] = params.cb
        t.root.icb_energy[:] = params_icb
        t.root.hb_energy[:] = params.hb
        t.root.ihb_energy[:] = params_ihb


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
        channel_membrane_potential=membrane_file,
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


def compute_divergence(config_base: str, pos: np.ndarray, mode: int = 0) -> Update:
    del mode
    with tb.open_file(config_base) as t:
        pot = t.root.input.potential
        if hasattr(pot, "cb_surf_membrane_potential"):
            cb_node_name = "cb_surf_membrane_potential"
            hb_node_name = "hb_surf_membrane_potential"
            cb_memb_shape = pot.cb_surf_membrane_potential.coeff.shape
            icb_memb_shape = pot.cb_surf_membrane_potential.coeff_inner.shape
            hb_memb_shape = pot.hb_surf_membrane_potential.coeff.shape
            ihb_memb_shape = pot.hb_surf_membrane_potential.coeff_inner.shape
        else:
            cb_node_name = "cb_membrane_potential"
            hb_node_name = "hb_membrane_potential"
            cb_memb_shape = pot.cb_membrane_potential.coeff.shape
            icb_memb_shape = (0,)
            hb_memb_shape = pot.hb_membrane_potential.coeff.shape
            ihb_memb_shape = (0,)

    cb_memb_size = int(np.prod(cb_memb_shape))
    icb_memb_size = int(np.prod(icb_memb_shape))
    hb_memb_size = int(np.prod(hb_memb_shape))
    ihb_memb_size = int(np.prod(ihb_memb_shape))

    engine = ue.Upside(config_base)
    contrast = Update([], [], [], [])

    for i in range(pos.shape[0]):
        engine.energy(pos[i])
        try:
            dp_cb_memb = engine.get_param_deriv((cb_memb_size + icb_memb_size,), cb_node_name)
            odp_cb_memb = dp_cb_memb[:cb_memb_size].reshape(cb_memb_shape)
            idp_cb_memb = (
                dp_cb_memb[cb_memb_size:].reshape(icb_memb_shape)
                if icb_memb_size
                else np.zeros(icb_memb_shape, dtype=np.float32)
            )
        except RuntimeError:
            dp_cb_memb = engine.get_param_deriv((cb_memb_size,), cb_node_name)
            odp_cb_memb = dp_cb_memb.reshape(cb_memb_shape)
            idp_cb_memb = np.zeros(icb_memb_shape, dtype=np.float32)

        try:
            dp_hb_memb = engine.get_param_deriv((hb_memb_size + ihb_memb_size,), hb_node_name)
            odp_hb_memb = dp_hb_memb[:hb_memb_size].reshape(hb_memb_shape)
            idp_hb_memb = (
                dp_hb_memb[hb_memb_size:].reshape(ihb_memb_shape)
                if ihb_memb_size
                else np.zeros(ihb_memb_shape, dtype=np.float32)
            )
        except RuntimeError:
            dp_hb_memb = engine.get_param_deriv((hb_memb_size,), hb_node_name)
            odp_hb_memb = dp_hb_memb.reshape(hb_memb_shape)
            idp_hb_memb = np.zeros(ihb_memb_shape, dtype=np.float32)

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
        return "local"
    return mode


def _run_worker_subprocess(
    state: dict,
    nm: str,
    direc: str,
    worker_argv: List[str],
) -> Tuple[sp.Popen, Optional[object]]:
    launch_mode = _choose_worker_launch(state["worker_launch"])
    n_threads = int(state["n_threads"])

    if launch_mode == "srun":
        cmd = [
            "srun",
            "--nodes=1",
            "--ntasks=1",
            f"--cpus-per-task={n_threads}",
            "--slurmd-debug=0",
            f"--output={direc}/{nm}.output_worker",
            state["worker_python"],
            *worker_argv,
        ]
        return sp.Popen(cmd, close_fds=True), None

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
                str(state["n_threads"]),
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

                rmsd[nm] = (divergence["rmsd_restrain"], divergence["rmsd"])
                com[nm] = (divergence["com_restrain"], divergence["com"])
                change.append(divergence["contrast"])
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
    }

    new_param_files = {k: os.path.join(direc, os.path.basename(x)) for k, x in initial_param_files.items()}
    new_param = param + step_update
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
    n_threads = int(sys.argv[14])
    project_root = sys.argv[15]
    ff_dir = sys.argv[16]
    del project_root

    frame_interval = max(1, int(sim_time / CONFIG.n_frame))
    path_code = f"{direc}/{code}"

    chain_break_compat = _prepare_chain_break_file(
        chain_break, os.path.join(direc, f"{code}.chain_break_compat.txt")
    )
    kwargs = _forcefield_kwargs(ff_dir, param_files["memb"], thickness, chain_break_compat)
    T = _temperature_schedule(n_threads)

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
    except Exception as exc:
        raise RuntimeError(f"CONFIG_FAIL: {exc}") from exc

    set1, set2 = gen_swap_set(1, len(T))
    swap_sets = [x for x in [set1, set2] if x]
    if not swap_sets:
        swap_sets = ru.swap_table2d(1, len(T))
    job = ru.run_upside(
        "",
        configs,
        sim_time,
        frame_interval,
        n_threads=n_threads,
        temperature=T,
        swap_sets=swap_sets,
        mc_interval=5.0,
        replica_interval=5.0,
        time_step=0.01,
        disable_z_recentering=True,
    )
    run_retcode = job.job.wait()
    if run_retcode != 0:
        # Upside occasionally exits non-zero after writing usable output.
        if not all(_has_output_frames(cfg) for cfg in configs[: min(2, len(configs))]):
            raise RuntimeError("RUN_FAIL")
        print(f"WARNING: run_upside returned {run_retcode}, continuing with existing output frames")

    divergence = {}
    equil_fraction = 0.25
    start = int(equil_fraction * CONFIG.n_frame)

    with tb.open_file(configs[0]) as t:
        pos_restrain = t.root.output.pos[start:, 0]
    with tb.open_file(configs[1]) as t:
        pos_free = t.root.output.pos[start:, 0]

    target = _load_native_positions(initial_structure_npy)
    divergence["rmsd_restrain"] = float(
        ru.traj_rmsd(pos_restrain[:, CONFIG.rmsd_k : -CONFIG.rmsd_k], target[CONFIG.rmsd_k : -CONFIG.rmsd_k]).mean()
    )
    divergence["rmsd"] = float(
        ru.traj_rmsd(pos_free[:, CONFIG.rmsd_k : -CONFIG.rmsd_k], target[CONFIG.rmsd_k : -CONFIG.rmsd_k]).mean()
    )
    divergence["com_restrain"] = float(pos_restrain[:, :, 2].mean())
    divergence["com"] = float(pos_free[:, :, 2].mean())

    all_pos = np.concatenate([pos_restrain, pos_free], axis=0)
    all_div = compute_divergence(configs[1], all_pos)
    n_pos0 = len(pos_restrain)
    divergence["contrast"] = Update(
        *[x[:n_pos0].mean(axis=0) - x[n_pos0:].mean(axis=0) for x in all_div]
    )
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

    if not __file__:
        raise RuntimeError("No __file__ available")
    state["worker_path"] = os.path.join(state["base_dir"], "ConDiv_mem.py")
    state["worker_python"] = CONFIG.worker_python
    shutil.copy(__file__, state["worker_path"])

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
        state["param"], state["init_param_files"] = get_init_param(state["init_dir"])
        with open(os.path.join(state["base_dir"], "condiv_init.pkl"), "wb") as fh:
            cp.dump((state["init_dir"], state["param"], state["init_param_files"]), fh, -1)
    else:
        with open(os.path.join(state["base_dir"], "condiv_init.pkl"), "rb") as fh:
            state["init_dir"], state["param"], state["init_param_files"] = cp.load(fh)

    state["initial_alpha"] = Update(0.02, 0.02, 0.02, 0.02) * CONFIG.alpha
    state["solver"] = AdamSolver(len(state["initial_alpha"]), alpha=state["initial_alpha"])
    state["sim_time"] = CONFIG.sim_time

    state["project_root"] = str(CONFIG.project_root)
    state["ff_dir"] = str(CONFIG.ff_dir)
    state["worker_launch"] = CONFIG.worker_launch
    state["n_threads"] = CONFIG.n_threads
    state["max_parallel_workers"] = CONFIG.max_parallel_workers

    print()
    print("Optimizing with solver", state["solver"])
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
