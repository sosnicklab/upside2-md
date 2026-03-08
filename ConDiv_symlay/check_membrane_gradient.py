#!/usr/bin/env python3
"""Finite-difference gradient checker for membrane ConDiv checkpoints."""

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
from typing import Tuple

import numpy as np
import tables as tb

DEFAULT_PROJECT_ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("CONDIV_PROJECT_ROOT", str(DEFAULT_PROJECT_ROOT))


def _load_module(worker_script: Path):
    spec = importlib.util.spec_from_file_location("ConDiv_mem", str(worker_script))
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load ConDiv_mem module from {worker_script}")
    module = importlib.util.module_from_spec(spec)
    sys.modules["ConDiv_mem"] = module
    spec.loader.exec_module(module)
    return module


def _guess_worker_script(checkpoint: Path) -> Path:
    candidates = [
        checkpoint.parent / "ConDiv_mem.py",
        checkpoint.parent.parent / "ConDiv_mem.py",
        checkpoint.parent.parent.parent / "ConDiv_mem.py",
    ]
    for cand in candidates:
        if cand.exists():
            return cand
    raise RuntimeError(
        f"Could not infer worker script for checkpoint {checkpoint}. "
        "Pass a checkpoint inside a ConDiv run directory."
    )


def _load_state(checkpoint: Path, worker_script: Path):
    _load_module(worker_script)
    with open(checkpoint, "rb") as fh:
        return cp.load(fh)


def _choose_target(state) -> Tuple[str, object]:
    candidates = []
    for minibatch in state["minibatches"]:
        for item in minibatch:
            candidates.append(item)
    if not candidates:
        raise RuntimeError("No targets available in checkpoint minibatches")
    return min(candidates, key=lambda x: x[1].n_res)


def _resolve_active_param_size(engine, node_name: str, coeff_size: int, inner_size: int) -> int:
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


def _select_fd_indices(analytic: np.ndarray, n_sample: int, rng: np.random.Generator) -> np.ndarray:
    if n_sample <= 0:
        return np.zeros((0,), dtype=np.int64)

    abs_an = np.abs(analytic.astype(np.float64, copy=False))
    finite_idx = np.where(np.isfinite(abs_an))[0]
    if finite_idx.size == 0:
        return np.zeros((0,), dtype=np.int64)

    if finite_idx.size <= n_sample:
        return finite_idx

    # Focus FD probes on high-signal coefficients to avoid near-zero noise dominance.
    pool_size = min(finite_idx.size, max(n_sample, n_sample * 8))
    if pool_size < finite_idx.size:
        top_local = np.argpartition(abs_an[finite_idx], -pool_size)[-pool_size:]
        candidate_idx = finite_idx[top_local]
    else:
        candidate_idx = finite_idx

    if candidate_idx.size <= n_sample:
        return candidate_idx
    return rng.choice(candidate_idx, size=n_sample, replace=False)


def _finite_diff_node(
    engine,
    pos: np.ndarray,
    node_name: str,
    active_size: int,
    eps: float,
    samples: int,
    rng: np.random.Generator,
):
    engine.energy(pos)
    analytic = engine.get_param_deriv((active_size,), node_name).astype(np.float64)
    base_param = engine.get_param((active_size,), node_name).astype(np.float64)

    n_sample = min(samples, active_size)
    idx = _select_fd_indices(analytic, n_sample, rng)

    fd_vals = []
    abs_err = []
    rel_err = []

    for i in idx:
        fd = np.nan
        an = float(analytic[i])
        eps_local = eps
        e_plus = np.nan
        e_minus = np.nan
        for _ in range(5):
            plus = base_param.copy()
            minus = base_param.copy()
            plus[i] += eps_local
            minus[i] -= eps_local

            engine.set_param(plus.astype(np.float32), node_name)
            e_plus = float(engine.energy(pos))
            engine.set_param(minus.astype(np.float32), node_name)
            e_minus = float(engine.energy(pos))

            if np.isfinite(e_plus) and np.isfinite(e_minus):
                fd = (e_plus - e_minus) / (2.0 * eps_local)
                break
            eps_local *= 0.1

        if not np.isfinite(fd):
            continue

        ae = abs(fd - an)
        re = ae / max(1e-6, abs(fd) + abs(an))

        fd_vals.append(fd)
        abs_err.append(ae)
        rel_err.append(re)

    engine.set_param(base_param.astype(np.float32), node_name)

    fd_vals = np.asarray(fd_vals, dtype=np.float64)
    abs_err = np.asarray(abs_err, dtype=np.float64)
    rel_err = np.asarray(rel_err, dtype=np.float64)
    n_used = int(fd_vals.size)

    return {
        "node": node_name,
        "active_size": int(active_size),
        "samples_requested": int(n_sample),
        "samples_used": n_used,
        "analytic_norm": float(np.linalg.norm(analytic)),
        "fd_norm": float(np.linalg.norm(fd_vals)) if n_used else float("nan"),
        "abs_err_mean": float(abs_err.mean()) if n_used else float("nan"),
        "abs_err_max": float(abs_err.max()) if n_used else float("nan"),
        "rel_err_median": float(np.median(rel_err)) if n_used else float("nan"),
        "rel_err_max": float(rel_err.max()) if n_used else float("nan"),
        "finite": bool(n_used > 0 and np.isfinite(analytic).all() and np.isfinite(fd_vals).all()),
    }


def _get_analytic_vector(engine, node_name: str, active_size: int):
    analytic = engine.get_param_deriv((active_size,), node_name).astype(np.float64)
    return analytic


def _perturb_membrane_file(
    membrane_path: Path,
    dataset_primary: str,
    dataset_secondary: str,
    primary_size: int,
    index: int,
    delta: float,
) -> None:
    with tb.open_file(str(membrane_path), "a") as t:
        if index < primary_size:
            arr = getattr(t.root, dataset_primary)[:]
            flat = arr.reshape(-1)
            flat[index] += delta
            getattr(t.root, dataset_primary)[:] = flat.reshape(arr.shape)
            return

        if not dataset_secondary:
            raise RuntimeError(
                f"Index {index} falls into secondary parameter block but no secondary dataset is available"
            )
        sec_idx = index - primary_size
        arr = getattr(t.root, dataset_secondary)[:]
        flat = arr.reshape(-1)
        if sec_idx < 0 or sec_idx >= flat.size:
            raise RuntimeError(
                f"Secondary index out of range: sec_idx={sec_idx}, size={flat.size}, index={index}"
            )
        flat[sec_idx] += delta
        getattr(t.root, dataset_secondary)[:] = flat.reshape(arr.shape)


def _finite_diff_by_rebuild(
    mod,
    target,
    ff_dir: str,
    chain_break_compat: str,
    init_npy: str,
    base_membrane: str,
    dataset_primary: str,
    dataset_secondary: str,
    primary_size: int,
    node_name: str,
    analytic: np.ndarray,
    active_size: int,
    eps: float,
    samples: int,
    rng: np.random.Generator,
    tmp_dir: Path,
):
    n_sample = min(samples, active_size)
    idx = _select_fd_indices(analytic, n_sample, rng)

    fd_vals = []
    abs_err = []
    rel_err = []

    for i in idx:
        e_vals = {}
        for sign, tag in ((+1.0, "plus"), (-1.0, "minus")):
            memb_path = tmp_dir / f"fd_{dataset_primary}_{i}_{tag}.memb.h5"
            shutil.copyfile(base_membrane, memb_path)
            _perturb_membrane_file(
                membrane_path=memb_path,
                dataset_primary=dataset_primary,
                dataset_secondary=dataset_secondary,
                primary_size=primary_size,
                index=int(i),
                delta=float(sign * eps),
            )

            cfg_path = tmp_dir / f"fd_{dataset_primary}_{i}_{tag}.cfg.h5"
            kwargs = mod._forcefield_kwargs(ff_dir, str(memb_path), str(target.thickness), chain_break_compat)
            mod.ru.upside_config(target.fasta, str(cfg_path), initial_structure=init_npy, **kwargs)
            mod.ru.advanced_config(
                str(cfg_path),
                restraint_groups=[f"0-{target.n_res - 1}"],
                restraint_spring_constant=mod.CONFIG.native_restraint_strength,
            )

            with tb.open_file(str(cfg_path)) as t:
                pos = t.root.input.pos[:, :, 0]
            eng = mod.ue.Upside(str(cfg_path))
            e_vals[tag] = float(eng.energy(pos))

        if not (np.isfinite(e_vals["plus"]) and np.isfinite(e_vals["minus"])):
            continue
        fd = (e_vals["plus"] - e_vals["minus"]) / (2.0 * eps)
        an = float(analytic[i])
        ae = abs(fd - an)
        re = ae / max(1e-6, abs(fd) + abs(an))

        fd_vals.append(fd)
        abs_err.append(ae)
        rel_err.append(re)

    fd_vals = np.asarray(fd_vals, dtype=np.float64)
    abs_err = np.asarray(abs_err, dtype=np.float64)
    rel_err = np.asarray(rel_err, dtype=np.float64)
    n_used = int(fd_vals.size)

    return {
        "node": node_name,
        "active_size": int(active_size),
        "samples_requested": int(n_sample),
        "samples_used": int(n_used),
        "analytic_norm": float(np.linalg.norm(analytic)),
        "fd_norm": float(np.linalg.norm(fd_vals)) if n_used else float("nan"),
        "abs_err_mean": float(abs_err.mean()) if n_used else float("nan"),
        "abs_err_max": float(abs_err.max()) if n_used else float("nan"),
        "rel_err_median": float(np.median(rel_err)) if n_used else float("nan"),
        "rel_err_max": float(rel_err.max()) if n_used else float("nan"),
        "finite": bool(n_used > 0 and np.isfinite(analytic).all() and np.isfinite(fd_vals).all()),
        "method": "config_rebuild_fd",
    }


def _sanity_from_state(state: dict):
    stats = state.get("last_gradient_stats", {})
    if not stats:
        return {
            "available": False,
            "pass": False,
            "reason": "last_gradient_stats missing from checkpoint",
        }

    numeric = []
    keys = []
    for k, v in stats.items():
        if isinstance(v, (int, float)):
            numeric.append(float(v))
            keys.append(k)

    finite = bool(np.isfinite(np.asarray(numeric)).all()) if numeric else False
    nontrivial = any(abs(v) > 1e-12 for v in numeric)
    return {
        "available": True,
        "pass": bool(finite and nontrivial),
        "finite": finite,
        "nontrivial": nontrivial,
        "keys": keys,
        "stats": stats,
    }


def main():
    parser = argparse.ArgumentParser(description="Check membrane gradient quality from a ConDiv checkpoint")
    parser.add_argument("--checkpoint", required=True, help="Checkpoint path (*.pkl)")
    parser.add_argument("--report", default="", help="JSON report output path")
    parser.add_argument("--fd-eps", type=float, default=1e-3)
    parser.add_argument("--fd-samples", type=int, default=24)
    parser.add_argument("--rel-median-threshold", type=float, default=1.2e-1)
    parser.add_argument("--rel-max-threshold", type=float, default=6.0e-1)
    parser.add_argument("--seed", type=int, default=1234)
    args = parser.parse_args()

    checkpoint = Path(args.checkpoint).expanduser().resolve()
    worker_script = _guess_worker_script(checkpoint)
    mod = _load_module(worker_script)
    state = _load_state(checkpoint, worker_script)

    base_dir = Path(state["base_dir"]).resolve()
    report_path = Path(args.report).expanduser().resolve() if args.report else (base_dir / "gradient_report.json")

    tmp_dir = base_dir / "gradient_check_tmp"
    if tmp_dir.exists():
        shutil.rmtree(tmp_dir)
    tmp_dir.mkdir(parents=True, exist_ok=True)

    _, target = _choose_target(state)

    temp_param = {"memb": str(tmp_dir / "membrane_fd.h5")}
    mod.expand_param(state["param"], state["init_param_files"], temp_param)

    ff_dir = state.get("ff_dir", str(mod.CONFIG.ff_dir))
    chain_break_compat = mod._prepare_chain_break_file(
        target.breakfile_path, str(tmp_dir / "gradient.chain_break_compat.txt")
    )
    kwargs = mod._forcefield_kwargs(ff_dir, temp_param["memb"], str(target.thickness), chain_break_compat)

    config_path = tmp_dir / "gradient_check.base.h5"
    init_npy = tmp_dir / "gradient_init.npy"
    np.save(str(init_npy), mod._load_native_positions(target.native_path))
    mod.ru.upside_config(target.fasta, str(config_path), initial_structure=str(init_npy), **kwargs)
    mod.ru.advanced_config(
        str(config_path),
        restraint_groups=[f"0-{target.n_res - 1}"],
        restraint_spring_constant=mod.CONFIG.native_restraint_strength,
    )

    with tb.open_file(str(config_path)) as t:
        pot = t.root.input.potential
        if hasattr(pot, "cb_surf_membrane_potential"):
            cb_node = "cb_surf_membrane_potential"
            hb_node = "hb_surf_membrane_potential"
            cb_shape = pot.cb_surf_membrane_potential.coeff.shape
            icb_shape = pot.cb_surf_membrane_potential.coeff_inner.shape
            hb_shape = pot.hb_surf_membrane_potential.coeff.shape
            ihb_shape = pot.hb_surf_membrane_potential.coeff_inner.shape
        else:
            cb_node = "cb_membrane_potential"
            hb_node = "hb_membrane_potential"
            cb_shape = pot.cb_membrane_potential.coeff.shape
            icb_shape = (0,)
            hb_shape = pot.hb_membrane_potential.coeff.shape
            ihb_shape = (0,)
        pos = t.root.input.pos[:, :, 0]

    cb_coeff = int(np.prod(cb_shape))
    cb_inner = int(np.prod(icb_shape))
    cb_total = cb_coeff + cb_inner
    hb_coeff = int(np.prod(hb_shape))
    hb_inner = int(np.prod(ihb_shape))
    hb_total = hb_coeff + hb_inner

    engine = mod.ue.Upside(str(config_path))
    rng = np.random.default_rng(args.seed)
    cb_active = _resolve_active_param_size(engine, cb_node, cb_coeff, cb_inner)
    hb_active = _resolve_active_param_size(engine, hb_node, hb_coeff, hb_inner)
    print(
        json.dumps(
            {
                "dimension_check": {
                    "cb": {"coeff": cb_coeff, "inner": cb_inner, "active": cb_active},
                    "hb": {"coeff": hb_coeff, "inner": hb_inner, "active": hb_active},
                }
            },
            indent=2,
        )
    )
    cb_analytic = _get_analytic_vector(engine, cb_node, cb_active)
    hb_analytic = _get_analytic_vector(engine, hb_node, hb_active)

    cb_report = _finite_diff_node(
        engine=engine,
        pos=pos,
        node_name=cb_node,
        active_size=cb_active,
        eps=args.fd_eps,
        samples=args.fd_samples,
        rng=rng,
    )
    if not cb_report["finite"]:
        cb_report = _finite_diff_by_rebuild(
            mod=mod,
            target=target,
            ff_dir=ff_dir,
            chain_break_compat=chain_break_compat,
            init_npy=str(init_npy),
            base_membrane=temp_param["memb"],
            dataset_primary="cb_energy",
            dataset_secondary="icb_energy" if cb_active == cb_total and cb_inner else "",
            primary_size=cb_coeff,
            node_name=cb_node,
            analytic=cb_analytic,
            active_size=cb_active,
            eps=args.fd_eps,
            samples=args.fd_samples,
            rng=rng,
            tmp_dir=tmp_dir,
        )

    hb_report = _finite_diff_node(
        engine=engine,
        pos=pos,
        node_name=hb_node,
        active_size=hb_active,
        eps=args.fd_eps,
        samples=args.fd_samples,
        rng=rng,
    )
    if not hb_report["finite"]:
        hb_report = _finite_diff_by_rebuild(
            mod=mod,
            target=target,
            ff_dir=ff_dir,
            chain_break_compat=chain_break_compat,
            init_npy=str(init_npy),
            base_membrane=temp_param["memb"],
            dataset_primary="hb_energy",
            dataset_secondary="ihb_energy" if hb_active == hb_total and hb_inner else "",
            primary_size=hb_coeff,
            node_name=hb_node,
            analytic=hb_analytic,
            active_size=hb_active,
            eps=args.fd_eps,
            samples=args.fd_samples,
            rng=rng,
            tmp_dir=tmp_dir,
        )

    sanity = _sanity_from_state(state)

    fd_pass = (
        cb_report["finite"]
        and hb_report["finite"]
        and cb_report["rel_err_median"] <= args.rel_median_threshold
        and hb_report["rel_err_median"] <= args.rel_median_threshold
        and cb_report["rel_err_max"] <= args.rel_max_threshold
        and hb_report["rel_err_max"] <= args.rel_max_threshold
    )

    report = {
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "checkpoint": str(checkpoint),
        "worker_script": str(worker_script),
        "target": {
            "fasta": target.fasta,
            "native_path": target.native_path,
            "thickness": target.thickness,
            "n_res": int(target.n_res),
        },
        "params": {
            "fd_eps": args.fd_eps,
            "fd_samples": args.fd_samples,
            "rel_median_threshold": args.rel_median_threshold,
            "rel_max_threshold": args.rel_max_threshold,
        },
        "dimension_check": {
            "cb": {"coeff": cb_coeff, "inner": cb_inner, "active": cb_active},
            "hb": {"coeff": hb_coeff, "inner": hb_inner, "active": hb_active},
        },
        "cb": cb_report,
        "hb": hb_report,
        "sanity": sanity,
        "sanity_pass": bool(sanity["pass"]),
        "fd_pass": bool(fd_pass),
        "pass": bool(fd_pass),
    }

    report_path.parent.mkdir(parents=True, exist_ok=True)
    with open(report_path, "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, sort_keys=True)

    print(json.dumps({"report": str(report_path), "pass": report["pass"]}, indent=2))
    if not report["pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
