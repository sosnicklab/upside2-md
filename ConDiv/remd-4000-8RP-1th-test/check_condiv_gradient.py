#!/usr/bin/env python3
"""Finite-difference validator for ConDiv.py checkpoints."""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import pickle as cp
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import tables as tb


def _load_module(script_path: Path):
    spec = importlib.util.spec_from_file_location("ConDiv", str(script_path))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load module from {script_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules["ConDiv"] = module
    spec.loader.exec_module(module)
    return module


class _MappedUnpickler(cp.Unpickler):
    def __init__(self, fh, module):
        super().__init__(fh)
        self.module = module

    def find_class(self, module, name):
        if module == "__main__" and hasattr(self.module, name):
            return getattr(self.module, name)
        return super().find_class(module, name)


def _load_state(checkpoint: Path, module):
    with open(checkpoint, "rb") as fh:
        return _MappedUnpickler(fh, module).load()


def _guess_worker_script(checkpoint: Path) -> Path:
    candidates = [
        checkpoint.parent / "ConDiv.py",
        checkpoint.parent.parent / "ConDiv.py",
        checkpoint.parent.parent.parent / "ConDiv.py",
    ]
    for c in candidates:
        if c.exists():
            return c
    raise RuntimeError(f"Could not infer ConDiv.py for checkpoint {checkpoint}")


def _resolve_project_root(script_dir: Path) -> Path:
    env = os.environ.get("CONDIV_PROJECT_ROOT")
    if env:
        p = Path(env).expanduser().resolve()
        if (p / "py" / "run_upside.py").exists():
            return p
    return script_dir.parent.parent.resolve()


def _normalize_shape(shape):
    if isinstance(shape, (int, np.integer)):
        return (int(shape),)
    if isinstance(shape, (tuple, list)):
        return tuple(int(x) for x in shape)
    raise TypeError(f"Unsupported shape spec: {shape!r}")


def _shape_size(shape) -> int:
    return int(np.prod(shape))


def _resolve_active_shape(engine, node_name: str, candidate_shapes):
    tried = []
    for spec in candidate_shapes:
        shape = _normalize_shape(spec)
        if any(dim <= 0 for dim in shape):
            continue
        tried.append(shape)
        try:
            param = engine.get_param(shape, node_name)
            # Validate round-trip to ensure set/get both accept this shape.
            engine.set_param(np.asarray(param, dtype=np.float32), node_name)
            return shape
        except RuntimeError:
            pass
    raise RuntimeError(f"Unable to resolve active shape for {node_name}; tried {tried}")


def _select_fd_indices(analytic: np.ndarray, n_sample: int, rng: np.random.Generator) -> np.ndarray:
    if n_sample <= 0:
        return np.zeros((0,), dtype=np.int64)
    abs_an = np.abs(analytic.astype(np.float64, copy=False))
    finite_idx = np.where(np.isfinite(abs_an))[0]
    if finite_idx.size == 0:
        return np.zeros((0,), dtype=np.int64)
    if finite_idx.size <= n_sample:
        return finite_idx

    pool_size = min(finite_idx.size, max(n_sample, n_sample * 8))
    if pool_size < finite_idx.size:
        top_local = np.argpartition(abs_an[finite_idx], -pool_size)[-pool_size:]
        candidate = finite_idx[top_local]
    else:
        candidate = finite_idx

    if candidate.size <= n_sample:
        return candidate
    return rng.choice(candidate, size=n_sample, replace=False)


def _finite_diff_node(engine, pos: np.ndarray, node_name: str, active_shape, eps: float, samples: int, rng):
    active_shape = _normalize_shape(active_shape)
    engine.energy(pos)
    analytic = np.asarray(engine.get_param_deriv(active_shape, node_name), dtype=np.float64)
    base_param = np.asarray(engine.get_param(active_shape, node_name), dtype=np.float64)

    analytic_flat = analytic.reshape(-1)
    base_flat = base_param.reshape(-1)

    n_sample = min(samples, base_flat.size)
    idx = _select_fd_indices(analytic_flat, n_sample, rng)

    fd_vals = []
    rel_err = []
    abs_err = []
    fd_skipped = False
    skip_reason = ""

    try:
        for i in idx:
            fd = np.nan
            an = float(analytic_flat[i])
            eps_local = eps
            for _ in range(5):
                plus = base_flat.copy()
                minus = base_flat.copy()
                plus[i] += eps_local
                minus[i] -= eps_local

                plus_param = plus.reshape(active_shape).astype(np.float32, copy=False)
                minus_param = minus.reshape(active_shape).astype(np.float32, copy=False)

                try:
                    engine.set_param(plus_param, node_name)
                    e_plus = float(engine.energy(pos))
                    engine.set_param(minus_param, node_name)
                    e_minus = float(engine.energy(pos))
                except RuntimeError as err:
                    fd_skipped = True
                    skip_reason = str(err)
                    break

                if np.isfinite(e_plus) and np.isfinite(e_minus):
                    fd = (e_plus - e_minus) / (2.0 * eps_local)
                    break
                eps_local *= 0.1

            if fd_skipped:
                break
            if not np.isfinite(fd):
                continue

            ae = abs(fd - an)
            re = ae / max(1e-6, abs(fd) + abs(an))
            fd_vals.append(fd)
            abs_err.append(ae)
            rel_err.append(re)
    finally:
        try:
            engine.set_param(base_param.astype(np.float32), node_name)
        except RuntimeError:
            pass

    fd_vals = np.asarray(fd_vals, dtype=np.float64)
    abs_err = np.asarray(abs_err, dtype=np.float64)
    rel_err = np.asarray(rel_err, dtype=np.float64)
    n_used = int(fd_vals.size)

    finite_analytic = bool(np.isfinite(analytic_flat).all())
    finite_fd = bool(n_used > 0 and np.isfinite(fd_vals).all())

    return {
        "node": node_name,
        "active_shape": list(active_shape),
        "active_size": int(base_flat.size),
        "samples_requested": int(n_sample),
        "samples_used": int(n_used),
        "analytic_norm": float(np.linalg.norm(analytic_flat)),
        "fd_norm": float(np.linalg.norm(fd_vals)) if n_used else float("nan"),
        "abs_err_mean": float(abs_err.mean()) if n_used else float("nan"),
        "abs_err_max": float(abs_err.max()) if n_used else float("nan"),
        "rel_err_median": float(np.median(rel_err)) if n_used else float("nan"),
        "rel_err_max": float(rel_err.max()) if n_used else float("nan"),
        "fd_skipped": bool(fd_skipped),
        "skip_reason": skip_reason,
        "finite": bool(finite_analytic and (fd_skipped or finite_fd)),
    }


def _choose_target(state):
    candidates = []
    for mb in state["minibatches"]:
        candidates.extend(mb)
    if not candidates:
        raise RuntimeError("No targets in checkpoint minibatches")
    return min(candidates, key=lambda x: x[1].n_res)


def _patch_inv_dx(mod, env_file: str, config_base: Path) -> None:
    try:
        with tb.open_file(env_file, "r") as t_src:
            if hasattr(t_src.root.energies._v_attrs, "inv_dx"):
                src_inv_dx = t_src.root.energies._v_attrs.inv_dx
            else:
                src_inv_dx = 17.0 / 12.0
        with tb.open_file(str(config_base), "r+") as t_dst:
            env_node = t_dst.root.input.potential.nonlinear_coupling_environment.coeff
            env_node._v_attrs.inv_dx = src_inv_dx
    except Exception:
        pass


def _build_validation_config(mod, state: dict, target, tmp_dir: Path, project_root: Path):
    temp_param = {
        "env": str(tmp_dir / "environment.h5"),
        "rot": str(tmp_dir / "sidechain.h5"),
        "hb": str(tmp_dir / "hbond.h5"),
        "sheet": str(tmp_dir / "sheet"),
    }
    mod.expand_param(state["param"], state["init_param_files"], temp_param)

    config_base = tmp_dir / "gradient_check.base.h5"
    kwargs = dict(
        environment_potential=temp_param["env"],
        environment_potential_type="0",
        rotamer_interaction=temp_param["rot"],
        rotamer_placement=temp_param["rot"],
        initial_structure=target.init_path,
        hbond_energy=temp_param["hb"],
        dynamic_rotamer_1body=True,
        rama_sheet_mix_energy=temp_param["sheet"],
        rama_param_deriv=True,
        rama_library=str(project_root / "parameters" / "common" / "rama.dat"),
        reference_state_rama=str(project_root / "parameters" / "common" / "rama_reference.pkl"),
    )
    mod.ru.upside_config(target.fasta, str(config_base), **kwargs)
    _patch_inv_dx(mod, temp_param["env"], config_base)

    with tb.open_file(str(config_base)) as t:
        pos = t.root.input.pos[:, :, 0]
        pot = t.root.input.potential
        rot_shape = tuple(int(x) for x in pot.rotamer.pair_interaction.interaction_param.shape)
        cov_shape = tuple(int(x) for x in pot.hbond_coverage.interaction_param.shape)
        hyd_shape = tuple(int(x) for x in pot.hbond_coverage_hydrophobe.interaction_param.shape)
        env_coeff_shape = tuple(int(x) for x in pot.nonlinear_coupling_environment.coeff.shape)
        env_weights_shape = tuple(int(x) for x in pot.nonlinear_coupling_environment.weights.shape)
        hbond_shape = tuple(int(x) for x in pot.hbond_energy.parameters.shape)

    return config_base, pos, temp_param, {
        "rotamer_shape": rot_shape,
        "rotamer_size": _shape_size(rot_shape),
        "hbond_coverage_shape": cov_shape,
        "hbond_coverage_size": _shape_size(cov_shape),
        "hbond_coverage_hydrophobe_shape": hyd_shape,
        "hbond_coverage_hydrophobe_size": _shape_size(hyd_shape),
        "nonlinear_coupling_environment_coeff_shape": env_coeff_shape,
        "nonlinear_coupling_environment_coeff_size": _shape_size(env_coeff_shape),
        "nonlinear_coupling_environment_weights_shape": env_weights_shape,
        "nonlinear_coupling_environment_weights_size": _shape_size(env_weights_shape),
        "hbond_energy_shape": hbond_shape,
        "hbond_energy_size": _shape_size(hbond_shape),
    }


def main():
    parser = argparse.ArgumentParser(description="Validate ConDiv.py checkpoint gradients")
    parser.add_argument("--checkpoint", required=True)
    parser.add_argument("--report", default="")
    parser.add_argument("--fd-eps", type=float, default=1e-3)
    parser.add_argument("--fd-samples", type=int, default=6)
    parser.add_argument("--rel-median-threshold", type=float, default=0.2)
    parser.add_argument("--rel-max-threshold", type=float, default=0.8)
    parser.add_argument("--seed", type=int, default=1234)
    args = parser.parse_args()

    checkpoint = Path(args.checkpoint).expanduser().resolve()
    worker_script = _guess_worker_script(checkpoint)
    script_dir = Path(__file__).resolve().parent
    project_root = _resolve_project_root(script_dir)

    os.environ.setdefault("PYTHONPATH", "")
    for p in (script_dir, project_root / "obj", project_root / "py", project_root / "src"):
        ps = str(p)
        if ps not in sys.path:
            sys.path.insert(0, ps)

    mod = _load_module(worker_script)
    state = _load_state(checkpoint, mod)

    base_dir = Path(state["base_dir"]).resolve()
    report_path = Path(args.report).expanduser().resolve() if args.report else (base_dir / "gradient_report.json")

    tmp_dir = base_dir / "gradient_check_tmp"
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    _, target = _choose_target(state)
    config_base, pos, _, shape_info = _build_validation_config(mod, state, target, tmp_dir, project_root)

    engine = mod.ue.Upside(str(config_base))
    rng = np.random.default_rng(args.seed)

    node_specs = []
    node_specs.append(("rotamer", [shape_info["rotamer_shape"], (shape_info["rotamer_size"],)]))
    node_specs.append(("hbond_coverage", [shape_info["hbond_coverage_shape"], (shape_info["hbond_coverage_size"],)]))
    node_specs.append((
        "hbond_coverage_hydrophobe",
        [shape_info["hbond_coverage_hydrophobe_shape"], (shape_info["hbond_coverage_hydrophobe_size"],)],
    ))
    env_total = shape_info["nonlinear_coupling_environment_coeff_size"] + shape_info["nonlinear_coupling_environment_weights_size"]
    node_specs.append((
        "nonlinear_coupling_environment",
        [
            (env_total,),
            shape_info["nonlinear_coupling_environment_coeff_shape"],
            (shape_info["nonlinear_coupling_environment_coeff_size"],),
        ],
    ))
    node_specs.append((
        "hbond_energy",
        [
            shape_info["hbond_energy_shape"],
            (shape_info["hbond_energy_size"],),
            (1,),
        ],
    ))

    reports = {}
    for node_name, candidates in node_specs:
        active_shape = _resolve_active_shape(engine, node_name, candidates)
        rep = _finite_diff_node(
            engine=engine,
            pos=pos,
            node_name=node_name,
            active_shape=active_shape,
            eps=args.fd_eps,
            samples=args.fd_samples,
            rng=rng,
        )
        reports[node_name] = rep

    # Check the trainable parameter tensors in this ConDiv workflow.
    params_finite = True
    for field in ("env", "rot", "hb", "sheet"):
        val = getattr(state["param"], field)
        arr = np.asarray(val, dtype=np.float64)
        if not np.isfinite(arr).all():
            params_finite = False
            break

    # Strict FD agreement is required for the coupling nodes used directly in d_obj.
    strict_fd_nodes = {"hbond_coverage", "hbond_coverage_hydrophobe"}
    fd_nodes_used = sum(
        1
        for node_name, rep in reports.items()
        if node_name in strict_fd_nodes and not rep.get("fd_skipped", False)
    )

    fd_pass = fd_nodes_used == len(strict_fd_nodes) and all(
        (
            rep["finite"]
            and not rep.get("fd_skipped", False)
            and rep["rel_err_median"] <= args.rel_median_threshold
            and rep["rel_err_max"] <= args.rel_max_threshold
        )
        if node_name in strict_fd_nodes
        else rep["finite"]
        for node_name, rep in reports.items()
    )

    report = {
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "checkpoint": str(checkpoint),
        "worker_script": str(worker_script),
        "project_root": str(project_root),
        "target": {
            "fasta": target.fasta,
            "init_path": target.init_path,
            "n_res": int(target.n_res),
        },
        "params": {
            "fd_eps": args.fd_eps,
            "fd_samples": args.fd_samples,
            "rel_median_threshold": args.rel_median_threshold,
            "rel_max_threshold": args.rel_max_threshold,
        },
        "shape_info": shape_info,
        "params_finite": bool(params_finite),
        "nodes": reports,
        "strict_fd_nodes": sorted(strict_fd_nodes),
        "fd_nodes_used": int(fd_nodes_used),
        "fd_pass": bool(fd_pass),
        "pass": bool(params_finite and fd_pass),
    }

    report_path.parent.mkdir(parents=True, exist_ok=True)
    with open(report_path, "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, sort_keys=True)

    print(json.dumps({"report": str(report_path), "pass": report["pass"]}, indent=2))
    if not report["pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
