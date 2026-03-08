#!/usr/bin/env python3
"""Progress logging and convergence checks for ConDiv_symlay training."""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


CONTINUE_EXIT_CODE = 10


@dataclass(frozen=True)
class Config:
    progress_log: Path
    status_json: Path
    grad_threshold: float
    update_threshold: float
    patience: int


def _load_json(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as fh:
        return json.load(fh)


def _read_progress_log(path: Path) -> List[Dict[str, Any]]:
    if not path.exists():
        return []

    records: List[Dict[str, Any]] = []
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line:
                records.append(json.loads(line))
    return records


def _append_progress_record(path: Path, record: Dict[str, Any]) -> Tuple[bool, List[Dict[str, Any]]]:
    records = _read_progress_log(path)
    if records and records[-1].get("checkpoint") == record["checkpoint"]:
        return False, records

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(record, sort_keys=True) + "\n")

    records.append(record)
    return True, records


def _parse_epoch_minibatch(checkpoint: Path) -> Tuple[Optional[int], Optional[int]]:
    match = re.match(r"epoch_(\d+)_minibatch_(\d+)$", checkpoint.parent.name)
    if not match:
        return None, None
    return int(match.group(1)), int(match.group(2))


def _combined_norm(stats: Dict[str, Any], prefix: str) -> float:
    total = 0.0
    for suffix in ("cb", "icb", "hb", "ihb"):
        value = float(stats.get(f"{prefix}_{suffix}", 0.0))
        total += value * value
    return math.sqrt(total)


def _build_record(checkpoint: Path, run_steps: int, slurm_job_id: str, resubmit_count: int) -> Dict[str, Any]:
    stats_path = checkpoint.parent / "gradient_stats.json"
    if not stats_path.exists():
        raise RuntimeError(f"Missing gradient stats file: {stats_path}")

    stats = _load_json(stats_path)
    epoch, minibatch = _parse_epoch_minibatch(checkpoint)
    record: Dict[str, Any] = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "checkpoint": str(checkpoint),
        "gradient_stats_path": str(stats_path),
        "epoch": epoch,
        "minibatch": minibatch,
        "run_steps": int(run_steps),
        "slurm_job_id": slurm_job_id,
        "resubmit_count": int(resubmit_count),
    }
    record.update(stats)
    record.setdefault("grad_norm_total", _combined_norm(record, "grad_norm"))
    record.setdefault("update_norm_total", _combined_norm(record, "update_norm"))
    return record


def _convergence_window(records: List[Dict[str, Any]], cfg: Config) -> Tuple[bool, List[Dict[str, Any]], str]:
    if len(records) < cfg.patience:
        return False, records[-cfg.patience :], "not_enough_history"

    window = records[-cfg.patience :]
    for rec in window:
        grad_norm = float(rec.get("grad_norm_total", float("inf")))
        update_norm = float(rec.get("update_norm_total", float("inf")))
        if grad_norm > cfg.grad_threshold or update_norm > cfg.update_threshold:
            return False, window, "threshold_not_met"

    return True, window, "converged"


def _status_payload(
    records: List[Dict[str, Any]],
    window: List[Dict[str, Any]],
    converged: bool,
    reason: str,
    cfg: Config,
) -> Dict[str, Any]:
    latest = records[-1] if records else {}
    return {
        "converged": bool(converged),
        "reason": reason,
        "n_records": int(len(records)),
        "patience": int(cfg.patience),
        "grad_threshold": float(cfg.grad_threshold),
        "update_threshold": float(cfg.update_threshold),
        "latest_checkpoint": latest.get("checkpoint"),
        "latest_grad_norm_total": latest.get("grad_norm_total"),
        "latest_update_norm_total": latest.get("update_norm_total"),
        "latest_epoch": latest.get("epoch"),
        "latest_minibatch": latest.get("minibatch"),
        "latest_slurm_job_id": latest.get("slurm_job_id"),
        "window": [
            {
                "checkpoint": rec.get("checkpoint"),
                "epoch": rec.get("epoch"),
                "minibatch": rec.get("minibatch"),
                "grad_norm_total": rec.get("grad_norm_total"),
                "update_norm_total": rec.get("update_norm_total"),
                "slurm_job_id": rec.get("slurm_job_id"),
            }
            for rec in window
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Append ConDiv training progress and check convergence")
    parser.add_argument("--base-dir", required=True, help="ConDiv run directory")
    parser.add_argument("--checkpoint", required=True, help="Latest checkpoint path")
    parser.add_argument("--progress-log", required=True, help="Progress log path (JSONL)")
    parser.add_argument("--status-json", required=True, help="Status JSON path")
    parser.add_argument("--grad-threshold", type=float, required=True, help="Convergence threshold for total gradient norm")
    parser.add_argument("--update-threshold", type=float, required=True, help="Convergence threshold for total update norm")
    parser.add_argument("--patience", type=int, required=True, help="Number of consecutive logged records required")
    parser.add_argument("--run-steps", type=int, required=True, help="Restart iterations executed by this job")
    parser.add_argument("--slurm-job-id", default="", help="Current Slurm job id")
    parser.add_argument("--resubmit-count", type=int, default=0, help="Auto-resubmission count for this job")
    args = parser.parse_args()

    cfg = Config(
        progress_log=Path(args.progress_log).expanduser().resolve(),
        status_json=Path(args.status_json).expanduser().resolve(),
        grad_threshold=float(args.grad_threshold),
        update_threshold=float(args.update_threshold),
        patience=int(args.patience),
    )

    checkpoint = Path(args.checkpoint).expanduser().resolve()
    record = _build_record(
        checkpoint=checkpoint,
        run_steps=int(args.run_steps),
        slurm_job_id=str(args.slurm_job_id),
        resubmit_count=int(args.resubmit_count),
    )
    appended, records = _append_progress_record(cfg.progress_log, record)
    converged, window, reason = _convergence_window(records, cfg)
    status = _status_payload(records, window, converged, reason, cfg)
    cfg.status_json.parent.mkdir(parents=True, exist_ok=True)
    with cfg.status_json.open("w", encoding="utf-8") as fh:
        json.dump(status, fh, indent=2, sort_keys=True)

    print("Training progress log:", cfg.progress_log)
    print("Training status json:", cfg.status_json)
    print("Progress record appended:", "yes" if appended else "no")
    print("Latest checkpoint:", status["latest_checkpoint"])
    print("Latest total grad norm:", status["latest_grad_norm_total"])
    print("Latest total update norm:", status["latest_update_norm_total"])
    print("Convergence status:", "converged" if converged else reason)

    return 0 if converged else CONTINUE_EXIT_CODE


if __name__ == "__main__":
    sys.exit(main())
