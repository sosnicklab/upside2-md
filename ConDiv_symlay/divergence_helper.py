#!/usr/bin/env python3
"""Minimal helper for ConDiv membrane divergence evaluation."""

from __future__ import annotations

import os
import pickle as cp
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Tuple

import numpy as np
import tables as tb


def _locate_project_root() -> Path:
    env_root = os.environ.get("CONDIV_PROJECT_ROOT")
    if env_root:
        root = Path(env_root).expanduser().resolve()
        if (root / "py" / "upside_engine.py").exists():
            return root

    starts = [Path(__file__).resolve(), Path.cwd().resolve()]
    seen: set[Path] = set()
    for start in starts:
        for candidate in [start.parent, *start.parents]:
            if candidate in seen:
                continue
            seen.add(candidate)
            if (candidate / "py" / "upside_engine.py").exists() and (candidate / "obj").exists():
                return candidate

    raise RuntimeError(
        "Could not locate project root containing py/upside_engine.py. "
        "Set CONDIV_PROJECT_ROOT to the repository root."
    )


PROJECT_ROOT = _locate_project_root()
for extra in (PROJECT_ROOT / "py", PROJECT_ROOT / "obj"):
    path_str = str(extra)
    if path_str not in sys.path:
        sys.path.insert(0, path_str)

import upside_engine as ue


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


def compute_divergence(config_base: str, pos: np.ndarray):
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
    contrast = [[], [], [], []]
    for frame in np.asarray(pos, dtype=np.float32):
        engine.energy(frame)

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

        contrast[0].append(odp_cb_memb)
        contrast[1].append(idp_cb_memb)
        contrast[2].append(odp_hb_memb)
        contrast[3].append(idp_hb_memb)

    return [np.asarray(term, dtype=np.float32) for term in contrast]


def main(argv: list[str]) -> int:
    if len(argv) != 4:
        print("Usage: divergence_helper.py <config.h5> <pos.pkl> <output.pkl>", file=sys.stderr)
        return 2

    os.environ["CONDIV_SKIP_UPSIDE_ENGINE_FREE"] = "1"
    with open(argv[2], "rb") as fh:
        pos = cp.load(fh)
    divergence = compute_divergence(argv[1], np.asarray(pos, dtype=np.float32))
    with open(argv[3], "wb") as fh:
        cp.dump(divergence, fh, protocol=4)
        fh.flush()
        os.fsync(fh.fileno())
    sys.stdout.flush()
    sys.stderr.flush()
    os._exit(0)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
