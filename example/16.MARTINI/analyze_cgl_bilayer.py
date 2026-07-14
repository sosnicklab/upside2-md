#!/usr/bin/env python3

import argparse
from pathlib import Path

import h5py
import numpy as np


COARSE_GRAIN_DIFFUSION_FACTOR = 14 * 4


def linear_fit(x, y):
    coeff = np.polyfit(x, y, 1)
    predicted = np.polyval(coeff, x)
    residual = np.sum((y - predicted) ** 2)
    total = np.sum((y - np.mean(y)) ** 2)
    r_squared = 1.0 - residual / total if total > 0.0 else 1.0
    return float(coeff[0]), float(coeff[1]), float(r_squared)


def lateral_diffusion(time, pos_xy, fit_min_time, fit_max_time):
    centered = pos_xy - np.mean(pos_xy, axis=1, keepdims=True)
    lag_time = time - time[0]
    lags = np.where(
        (lag_time >= float(fit_min_time))
        & (lag_time <= float(fit_max_time))
    )[0]
    if lags.size < 3:
        raise ValueError(
            f"trajectory does not contain three lags in the "
            f"{fit_min_time:g}-{fit_max_time:g} tu fit window"
        )
    msd = np.empty(lags.size, dtype=np.float64)
    for i, lag in enumerate(lags):
        delta = centered[lag:] - centered[:-lag]
        msd[i] = np.mean(np.sum(delta * delta, axis=2))

    fit_time = lag_time[lags]
    slope, intercept, r_squared = linear_fit(fit_time, msd)
    anomalous_exponent, _, _ = linear_fit(np.log(fit_time), np.log(msd))
    return slope / 4.0, r_squared, anomalous_exponent


def single_origin_diffusion(time, pos_xy, origin_fraction):
    centered = pos_xy - np.mean(pos_xy, axis=1, keepdims=True)
    origin = int(round(float(origin_fraction) * (time.size - 1)))
    elapsed = float(time[-1] - time[origin])
    if elapsed <= 0.0:
        raise ValueError("single-origin diffusion requires positive elapsed time")
    displacement = centered[-1] - centered[origin]
    return float(np.mean(np.sum(displacement * displacement, axis=1)) / (4.0 * elapsed))


def unwrap_positions(pos, box_lengths):
    if box_lengths is None:
        return np.asarray(pos, dtype=np.float64).copy()
    unwrapped = np.asarray(pos, dtype=np.float64).copy()
    box = np.asarray(box_lengths, dtype=np.float64)
    if box.ndim == 1:
        box = np.repeat(box[None, :], unwrapped.shape[0], axis=0)
    if box.shape != (unwrapped.shape[0], 3):
        raise ValueError("box lengths must have shape 3 or n_frame x 3")
    delta = unwrapped[1:] - unwrapped[:-1]
    step_box = box[1:, None, :]
    delta -= step_box * np.round(delta / step_box)
    unwrapped[1:] = unwrapped[0] + np.cumsum(delta, axis=0)
    return unwrapped


def lateral_uniformity(pos, box_lengths, leaflet_mask, grid_size=6):
    modes = np.asarray(((1, 0), (0, 1), (1, 1), (1, -1)), dtype=np.float64)
    structure = []
    occupancy_cv = []
    empty_fraction = []
    nearest = []
    for frame, box in zip(pos, box_lengths):
        xy = np.mod(frame[leaflet_mask, :2], box[:2])
        phase = 2.0 * np.pi * (
            xy[:, 0, None] * modes[None, :, 0] / box[0]
            + xy[:, 1, None] * modes[None, :, 1] / box[1]
        )
        structure.extend(np.abs(np.mean(np.exp(1j * phase), axis=0)))

        cell = np.floor(xy / box[:2] * int(grid_size)).astype(int)
        cell = np.clip(cell, 0, int(grid_size) - 1)
        counts = np.zeros((int(grid_size), int(grid_size)), dtype=np.float64)
        np.add.at(counts, (cell[:, 0], cell[:, 1]), 1.0)
        occupancy_cv.append(float(np.std(counts) / np.mean(counts)))
        empty_fraction.append(float(np.mean(counts == 0.0)))

        delta = xy[:, None, :] - xy[None, :, :]
        delta -= box[None, None, :2] * np.round(
            delta / box[None, None, :2]
        )
        distance = np.sqrt(np.sum(delta * delta, axis=2))
        np.fill_diagonal(distance, np.inf)
        nearest.extend(np.min(distance, axis=1))
    return {
        "structure": np.asarray(structure, dtype=np.float64),
        "occupancy_cv": np.asarray(occupancy_cv, dtype=np.float64),
        "empty_fraction": np.asarray(empty_fraction, dtype=np.float64),
        "nearest": np.asarray(nearest, dtype=np.float64),
    }


def analyze(path, fit_min_time, fit_max_time):
    with h5py.File(path, "r") as h5:
        time = np.asarray(h5["/output/time"][:], dtype=np.float64)
        pos = np.asarray(h5["/output/pos"][:, 0], dtype=np.float64)
        if "/input/potential/compose_vector6d/elem_index" in h5:
            cgl_indices = np.asarray(
                h5["/input/potential/compose_vector6d/elem_index"][:], dtype=int
            )
            pos = pos[:, cgl_indices]
        orientation = np.asarray(h5["/output/cgl_orientation"][:, 0], dtype=np.float64)
        compaction = (
            np.asarray(h5["/output/cgl_compaction"][:, 0], dtype=np.float64)
            if "/output/cgl_compaction" in h5
            else None
        )
        q_group = h5.get("/input/potential/cgl_response_state")
        if q_group is None:
            q_group = h5.get("/input/potential/cgl_compaction_state")
        q_min = float(q_group.attrs["self_coord_min_ang"]) if q_group else None
        q_max = float(q_group.attrs["self_coord_max_ang"]) if q_group else None
        box_lengths = None
        for group_path in (
            "/input/potential/cg_lipid_pair",
            "/input/potential/martini_potential",
        ):
            if group_path not in h5:
                continue
            attrs = h5[group_path].attrs
            if all(name in attrs for name in ("x_len", "y_len", "z_len")):
                candidate = np.asarray(
                    [attrs["x_len"], attrs["y_len"], attrs["z_len"]],
                    dtype=np.float64,
                )
                if np.all(candidate > 0.0):
                    box_lengths = candidate
                    break
        if "/output/box" in h5:
            output_box = np.asarray(h5["/output/box"][:], dtype=np.float64)
            if output_box.shape == (time.size, 3) and np.all(output_box > 0.0):
                box_lengths = output_box

    if not (
        np.all(np.isfinite(time))
        and np.all(np.isfinite(pos))
        and np.all(np.isfinite(orientation))
        and (compaction is None or np.all(np.isfinite(compaction)))
    ):
        raise RuntimeError(f"{path}: non-finite trajectory values")

    if box_lengths is None:
        raise ValueError(f"{path}: no valid periodic box lengths")
    if np.asarray(box_lengths).ndim == 1:
        box_lengths = np.repeat(
            np.asarray(box_lengths, dtype=np.float64)[None, :],
            time.size,
            axis=0,
        )
    unwrapped_pos = unwrap_positions(pos, box_lengths)
    full_time = time
    full_unwrapped_pos = unwrapped_pos
    first = int(np.floor(0.20 * time.size))
    time = time[first:]
    pos = pos[first:]
    box_lengths = box_lengths[first:]
    unwrapped_pos = unwrapped_pos[first:]
    orientation = orientation[first:]
    if compaction is not None:
        compaction = compaction[first:]

    initial_z = pos[0, :, 2]
    middle_z = float(np.median(initial_z))
    upper = initial_z >= middle_z
    lower = ~upper
    current_upper = pos[:, :, 2] >= np.median(pos[:, :, 2], axis=1)[:, None]
    flip_fraction = np.mean(current_upper != upper[None, :], axis=1)
    separation = np.mean(pos[:, upper, 2], axis=1) - np.mean(pos[:, lower, 2], axis=1)
    separation_slope, _, separation_r2 = linear_fit(time - time[0], separation)

    norm = np.linalg.norm(orientation, axis=2)
    abs_nz = np.abs(orientation[:, :, 2]) / np.maximum(norm, 1.0e-12)
    mean_abs_nz = np.mean(abs_nz, axis=1)
    order_slope, _, order_r2 = linear_fit(time - time[0], mean_abs_nz)

    diffusion, diffusion_r2, anomalous_exponent = lateral_diffusion(
        time,
        unwrapped_pos[:, :, :2],
        fit_min_time,
        fit_max_time,
    )
    single_origin = {
        label: single_origin_diffusion(full_time, full_unwrapped_pos[:, :, :2], fraction)
        for label, fraction in (
            ("full", 0.0),
            ("after_quarter", 0.25),
            ("late_half", 0.5),
            ("last_quarter", 0.75),
        )
    }
    uniformity = [
        lateral_uniformity(pos, box_lengths, mask)
        for mask in (upper, lower)
    ]
    print(path)
    print(f"  frames={time.size} retained_time={time[-1] - time[0]:.6g} tu")
    print(
        f"  D_multi_origin_raw={diffusion:.6g} A^2/T_up "
        f"D_x56={COARSE_GRAIN_DIFFUSION_FACTOR * diffusion:.6g} um^2/s "
        f"R2={diffusion_r2:.6f} alpha={anomalous_exponent:.6g} "
        f"fit=[{fit_min_time:g},{fit_max_time:g}] tu"
    )
    print(
        "  D_single_origin_raw "
        + " ".join(f"{label}={value:.6g}" for label, value in single_origin.items())
        + " A^2/T_up"
    )
    print(
        "  D_single_origin_x56 "
        + " ".join(
            f"{label}={COARSE_GRAIN_DIFFUSION_FACTOR * value:.6g}"
            for label, value in single_origin.items()
        )
        + " um^2/s"
    )
    print(
        f"  leaflet_separation mean={np.mean(separation):.6g} "
        f"initial={separation[0]:.6g} final={separation[-1]:.6g} "
        f"slope={separation_slope:.6g} A/tu R2={separation_r2:.6f}"
    )
    print(
        f"  abs_nz mean={np.mean(abs_nz):.6g} p05={np.quantile(abs_nz, 0.05):.6g} "
        f"final_mean={mean_abs_nz[-1]:.6g} slope={order_slope:.6g}/tu "
        f"R2={order_r2:.6f}"
    )
    for label, metrics in zip(("upper", "lower"), uniformity):
        print(
            f"  uniformity {label} low_k_mean/max="
            f"{np.mean(metrics['structure']):.6g}/{np.max(metrics['structure']):.6g} "
            f"occupancy6_cv_mean/max="
            f"{np.mean(metrics['occupancy_cv']):.6g}/"
            f"{np.max(metrics['occupancy_cv']):.6g} "
            f"empty6_mean/max={np.mean(metrics['empty_fraction']):.6g}/"
            f"{np.max(metrics['empty_fraction']):.6g}"
        )
        print(
            f"  nearest_xy {label} min/p05/mean/cv="
            f"{np.min(metrics['nearest']):.6g}/"
            f"{np.quantile(metrics['nearest'], 0.05):.6g}/"
            f"{np.mean(metrics['nearest']):.6g}/"
            f"{np.std(metrics['nearest']) / np.mean(metrics['nearest']):.6g} A"
        )
    print(
        f"  leaflet_flip_fraction mean/max/final={np.mean(flip_fraction):.6g}/"
        f"{np.max(flip_fraction):.6g}/{flip_fraction[-1]:.6g}"
    )
    if compaction is not None:
        edge_width = 0.01 * (q_max - q_min)
        print(
            f"  compaction mean={np.mean(compaction):.6g} std={np.std(compaction):.6g} "
            f"range=[{np.min(compaction):.6g},{np.max(compaction):.6g}] "
            f"edge_fraction={np.mean((compaction <= q_min + edge_width) | (compaction >= q_max - edge_width)):.6g}"
        )
    else:
        print("  compaction=not_present")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("trajectory", nargs="+", type=Path)
    parser.add_argument("--fit-min-time", type=float, default=1.0)
    parser.add_argument("--fit-max-time", type=float, default=5.0)
    args = parser.parse_args()
    if args.fit_min_time <= 0.0 or args.fit_max_time <= args.fit_min_time:
        parser.error("diffusion fit requires 0 < fit-min-time < fit-max-time")
    for trajectory in args.trajectory:
        analyze(trajectory, args.fit_min_time, args.fit_max_time)


if __name__ == "__main__":
    main()
