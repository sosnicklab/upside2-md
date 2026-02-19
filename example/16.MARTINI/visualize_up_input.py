#!../../.venv/bin/python
import argparse
from dataclasses import dataclass
import os
import tempfile
import numpy as np
import h5py

if "MPLCONFIGDIR" not in os.environ:
    os.environ["MPLCONFIGDIR"] = os.path.join(tempfile.gettempdir(), "matplotlib")
if "XDG_CACHE_HOME" not in os.environ:
    os.environ["XDG_CACHE_HOME"] = os.path.join(tempfile.gettempdir(), "xdg-cache")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)
os.makedirs(os.path.join(os.environ["XDG_CACHE_HOME"], "fontconfig"), exist_ok=True)
if "DISPLAY" not in os.environ and "MPLBACKEND" not in os.environ:
    os.environ["MPLBACKEND"] = "Agg"

import matplotlib.pyplot as plt


@dataclass
class Config:
    up_file: str
    output: str | None
    color_by: str
    point_size: float


def _default_output_path(up_file: str) -> str:
    base, _ = os.path.splitext(os.path.basename(up_file))
    out_name = f"{base}_input_viz.png"
    return os.path.join(os.path.dirname(up_file) or ".", out_name)


def _decode_str_array(arr: np.ndarray) -> np.ndarray:
    if arr.dtype.kind == "S":
        return np.array([x.decode("utf-8", errors="ignore").strip() for x in arr], dtype=object)
    return arr.astype(object)


def _load_up(path: str):
    with h5py.File(path, "r") as f:
        if "/input/pos" not in f:
            raise ValueError("Missing dataset: /input/pos")

        pos = f["/input/pos"][:]
        if pos.ndim == 3 and pos.shape[2] == 1:
            pos = pos[:, :, 0]
        elif pos.ndim != 2 or pos.shape[1] != 3:
            raise ValueError(f"Unexpected /input/pos shape: {pos.shape}")

        n = pos.shape[0]

        mass = f["/input/mass"][:] if "/input/mass" in f else np.ones(n, dtype=float)
        charges = f["/input/charges"][:] if "/input/charges" in f else np.zeros(n, dtype=float)
        residue_ids = (
            f["/input/residue_ids"][:] if "/input/residue_ids" in f else np.arange(n, dtype=int)
        )

        atom_type = None
        if "/input/type" in f:
            atom_type = _decode_str_array(f["/input/type"][:])

    return pos, mass, charges, residue_ids, atom_type


def _build_colors(cfg: Config, mass, charges, residue_ids, atom_type):
    if cfg.color_by == "charge":
        values = charges.astype(float)
        cmap = "coolwarm"
        label = "charge"
        return values, cmap, label

    if cfg.color_by == "mass":
        values = mass.astype(float)
        cmap = "viridis"
        label = "mass"
        return values, cmap, label

    if cfg.color_by == "residue":
        values = residue_ids.astype(float)
        cmap = "turbo"
        label = "residue_id"
        return values, cmap, label

    if cfg.color_by == "type":
        if atom_type is None:
            raise ValueError("Requested --color-by type, but /input/type is missing.")
        uniq = {k: i for i, k in enumerate(np.unique(atom_type))}
        values = np.array([uniq[x] for x in atom_type], dtype=float)
        cmap = "tab20"
        label = "type_index"
        return values, cmap, label

    raise ValueError(f"Unknown color mode: {cfg.color_by}")


def main():
    p = argparse.ArgumentParser(description="Visualize UPSIDE .up /input datasets")
    p.add_argument(
        "up_file",
        nargs="?",
        default="inputs/1rkl.up",
        help="Path to .up file (default: inputs/1rkl.up)",
    )
    p.add_argument(
        "--output",
        default=None,
        help="Output PNG path. If omitted, saves to <input_basename>_input_viz.png",
    )
    p.add_argument(
        "--color-by",
        choices=["charge", "mass", "residue", "type"],
        default="charge",
        help="Color scatter points by this field",
    )
    p.add_argument("--point-size", type=float, default=8.0)
    args = p.parse_args()

    cfg = Config(
        up_file=args.up_file,
        output=args.output,
        color_by=args.color_by,
        point_size=args.point_size,
    )

    pos, mass, charges, residue_ids, atom_type = _load_up(cfg.up_file)
    color_values, cmap, color_label = _build_colors(cfg, mass, charges, residue_ids, atom_type)

    print(f"Loaded: {cfg.up_file}")
    print(f"Atoms: {pos.shape[0]}")
    print(f"mass: shape={mass.shape}, min={mass.min():.4g}, max={mass.max():.4g}")
    print(
        f"charges: shape={charges.shape}, nonzero={(np.abs(charges) > 0).sum()}, "
        f"sum={charges.sum():.4g}"
    )
    print(f"pos: shape={pos.shape}, x[{pos[:,0].min():.3f},{pos[:,0].max():.3f}] "
          f"y[{pos[:,1].min():.3f},{pos[:,1].max():.3f}] z[{pos[:,2].min():.3f},{pos[:,2].max():.3f}]")

    fig, axes = plt.subplots(2, 2, figsize=(12, 10), constrained_layout=True)

    sc1 = axes[0, 0].scatter(
        pos[:, 0], pos[:, 1], c=color_values, s=cfg.point_size, cmap=cmap, alpha=0.8
    )
    axes[0, 0].set_title("XY projection")
    axes[0, 0].set_xlabel("X")
    axes[0, 0].set_ylabel("Y")
    axes[0, 0].set_aspect("equal", adjustable="box")

    sc2 = axes[0, 1].scatter(
        pos[:, 0], pos[:, 2], c=color_values, s=cfg.point_size, cmap=cmap, alpha=0.8
    )
    axes[0, 1].set_title("XZ projection")
    axes[0, 1].set_xlabel("X")
    axes[0, 1].set_ylabel("Z")
    axes[0, 1].set_aspect("equal", adjustable="box")

    axes[1, 0].hist(charges, bins=60, color="steelblue", edgecolor="black")
    axes[1, 0].set_title("Charge distribution")
    axes[1, 0].set_xlabel("charge")
    axes[1, 0].set_ylabel("count")

    axes[1, 1].hist(mass, bins=60, color="darkorange", edgecolor="black")
    axes[1, 1].set_title("Mass distribution")
    axes[1, 1].set_xlabel("mass")
    axes[1, 1].set_ylabel("count")

    cbar = fig.colorbar(sc2, ax=axes[:2, :2], shrink=0.8)
    cbar.set_label(color_label)

    output_path = cfg.output or _default_output_path(cfg.up_file)
    fig.savefig(output_path, dpi=200)
    print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()
