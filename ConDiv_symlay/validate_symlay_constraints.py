#!/usr/bin/env python3
"""Validate symmetry and multilayer constraints for ConDiv_symlay checkpoints or membrane files."""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import pickle as cp
import shutil
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("CONDIV_PROJECT_ROOT", str(Path(__file__).resolve().parents[1]))

from symlay_utils import (
    DEFAULT_DENSE_GRID_SIZE,
    DEFAULT_SUPPORT_MARGIN,
    analyze_membrane_constraints,
    load_layer_manifest,
    repo_relative,
)


def _load_module(worker_script: Path):
    script_dir = str(worker_script.parent)
    if script_dir not in sys.path:
        sys.path.insert(0, script_dir)
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
        "Pass a checkpoint inside a ConDiv_symlay run directory."
    )


def _materialize_from_checkpoint(checkpoint: Path, output_dir: Path) -> Path:
    worker_script = _guess_worker_script(checkpoint)
    mod = _load_module(worker_script)
    with checkpoint.open("rb") as fh:
        state = cp.load(fh)

    output_dir.mkdir(parents=True, exist_ok=True)
    temp_membrane = output_dir / "checkpoint_materialized_membrane.h5"
    mod.expand_param(state["param"], state["init_param_files"], {"memb": str(temp_membrane)})
    return temp_membrane


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Validate ConDiv_symlay symmetric multilayer constraints")
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--checkpoint", type=Path, help="ConDiv_symlay checkpoint to validate")
    source.add_argument("--membrane-h5", type=Path, help="Direct membrane.h5 file to validate")
    parser.add_argument(
        "--layer-manifest",
        type=Path,
        required=True,
        help="Layer manifest JSON produced by build_layer_manifest.py",
    )
    parser.add_argument("--report", type=Path, help="Optional JSON report output path")
    parser.add_argument(
        "--dense-grid-size",
        type=int,
        default=DEFAULT_DENSE_GRID_SIZE,
        help="Dense grid size used for projection-residual validation",
    )
    parser.add_argument(
        "--support-margin",
        type=float,
        default=DEFAULT_SUPPORT_MARGIN,
        help="Support margin used to compute the required symmetric support",
    )
    parser.add_argument(
        "--symmetry-threshold",
        type=float,
        default=1e-1,
        help="Maximum allowed symmetry residual",
    )
    parser.add_argument(
        "--projection-threshold",
        type=float,
        default=2.5e-1,
        help="Maximum allowed projection residual",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    manifest = load_layer_manifest(args.layer_manifest)
    os.environ["CONDIV_SYMLAY_LAYER_MANIFEST"] = str(args.layer_manifest.expanduser().resolve())

    temp_root = None
    if args.checkpoint is not None:
        temp_root = Path(tempfile.mkdtemp(prefix="condiv_symlay_validate_"))
        membrane_path = _materialize_from_checkpoint(args.checkpoint.expanduser().resolve(), temp_root)
    else:
        membrane_path = args.membrane_h5.expanduser().resolve()

    try:
        report = analyze_membrane_constraints(
            membrane_path,
            manifest,
            dense_grid_size=args.dense_grid_size,
            support_margin=args.support_margin,
        )
    finally:
        if temp_root is not None:
            shutil.rmtree(temp_root, ignore_errors=True)

    max_symmetry = max(
        float(report["projection_summary"][key]["max_symmetry_residual"]) for key in ("cb", "icb", "hb", "ihb")
    )
    max_projection = max(
        float(report["projection_summary"][key]["max_projection_residual"]) for key in ("cb", "icb", "hb", "ihb")
    )
    support_ok = all(bool(v) for v in report["support_checks"].values())
    passed = support_ok and max_symmetry <= args.symmetry_threshold and max_projection <= args.projection_threshold
    report["thresholds"] = {
        "symmetry_threshold": float(args.symmetry_threshold),
        "projection_threshold": float(args.projection_threshold),
    }
    report["pass"] = bool(passed)

    if args.report:
        report_path = args.report.expanduser().resolve()
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print("ConDiv_symlay constraint validation")
    print(f"  membrane: {repo_relative(membrane_path)}")
    print(f"  manifest: {repo_relative(args.layer_manifest)}")
    print(f"  pass: {report['pass']}")
    print(f"  max symmetry residual: {max_symmetry:.6g}")
    print(f"  max projection residual: {max_projection:.6g}")
    print(
        "  support checks: "
        + ", ".join(f"{key}={int(bool(val))}" for key, val in sorted(report["support_checks"].items()))
    )

    if not report["pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
