#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np


def ca_mask(atom_names: list[str]) -> np.ndarray:
    return np.asarray([name == "CA" for name in atom_names], dtype=bool)


def radius_of_gyration(coords: np.ndarray) -> np.ndarray:
    center = np.mean(coords, axis=1, keepdims=True)
    sq = np.sum((coords - center) ** 2, axis=2)
    return np.sqrt(np.mean(sq, axis=1))


def line_bend_rms(ca_coords: np.ndarray) -> np.ndarray:
    centered = ca_coords - np.mean(ca_coords, axis=1, keepdims=True)
    out = np.zeros(centered.shape[0], dtype=np.float64)
    for i, frame in enumerate(centered):
        _, _, vh = np.linalg.svd(frame, full_matrices=False)
        axis = vh[0]
        proj = np.outer(frame @ axis, axis)
        resid = frame - proj
        out[i] = np.sqrt(np.mean(np.sum(resid**2, axis=1)))
    return out


def summarize(values: np.ndarray) -> dict[str, float]:
    arr = np.asarray(values, dtype=np.float64)
    return {
        "mean": float(np.mean(arr)),
        "p05": float(np.quantile(arr, 0.05)),
        "p50": float(np.quantile(arr, 0.50)),
        "p95": float(np.quantile(arr, 0.95)),
        "min": float(np.min(arr)),
        "max": float(np.max(arr)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--up", required=True)
    parser.add_argument("--start-frac", type=float, default=0.5)
    args = parser.parse_args()

    up = Path(args.up).expanduser().resolve()
    with h5py.File(up, "r") as h5:
        atom_names = [
            x.decode("utf-8", errors="ignore") if isinstance(x, bytes) else str(x)
            for x in h5["/input/atom_names"][:]
        ]
        pos = np.asarray(h5["/output/pos"][:, 0, :, :], dtype=np.float64)
        start = min(max(int(pos.shape[0] * float(args.start_frac)), 0), pos.shape[0] - 1)
        pos = pos[start:]
        ca = ca_mask(atom_names)
        ca_pos = pos[:, ca, :]
        hbonds = np.sum(np.asarray(h5["/output/hbond"][start:], dtype=np.float64), axis=1)
        q = None
        if "/output/cgl_compaction" in h5:
            q = np.asarray(h5["/output/cgl_compaction"][start:, 0, :], dtype=np.float64)

    result = {
        "up": str(up),
        "frames": int(pos.shape[0]),
        "late_start_frame": int(start),
        "hbonds": summarize(hbonds),
        "rg_ca": summarize(radius_of_gyration(ca_pos)),
        "bend_rms_ca": summarize(line_bend_rms(ca_pos)),
        "end_to_end_ca": summarize(np.linalg.norm(ca_pos[:, -1, :] - ca_pos[:, 0, :], axis=1)),
    }
    if q is not None:
        flat_q = q.reshape(-1)
        result["cgl_compaction_q"] = summarize(flat_q)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
