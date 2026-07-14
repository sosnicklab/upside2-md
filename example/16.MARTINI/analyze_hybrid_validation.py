#!/usr/bin/env python3

import argparse
import shlex
from pathlib import Path

import h5py
import mdtraj as md
import numpy as np

from analyze_cgl_bilayer import lateral_uniformity


def decode(values):
    return np.asarray(
        [v.decode("ascii") if isinstance(v, bytes) else str(v) for v in values]
    )


def integration_stage_time_step(h5):
    time = h5["output/time"]
    if "integration_stage_time_step" in time.attrs:
        return float(time.attrs["integration_stage_time_step"])
    invocation = h5["output"].attrs.get("invocation", b"")
    if isinstance(invocation, bytes):
        invocation = invocation.decode("ascii")
    arguments = shlex.split(str(invocation))
    try:
        return float(arguments[arguments.index("--time-step") + 1])
    except (ValueError, IndexError) as exc:
        raise ValueError(
            "trajectory does not record its integration-stage timestep"
        ) from exc


def linear_fit(x, y):
    slope, intercept = np.polyfit(x, y, 1)
    predicted = slope * x + intercept
    total = np.sum((y - np.mean(y)) ** 2)
    residual = np.sum((y - predicted) ** 2)
    r_squared = 1.0 - residual / total if total > 0.0 else 1.0
    return float(slope), float(r_squared)


def kabsch_rotation(coords, reference):
    centered = coords - np.mean(coords, axis=0)
    ref_centered = reference - np.mean(reference, axis=0)
    u, _, vt = np.linalg.svd(centered.T @ ref_centered)
    handedness = np.linalg.det(u @ vt)
    return u @ np.diag([1.0, 1.0, handedness]) @ vt


def kabsch_rmsd(coords, reference):
    centered = coords - np.mean(coords, axis=0)
    ref_centered = reference - np.mean(reference, axis=0)
    rotation = kabsch_rotation(coords, reference)
    delta = centered @ rotation - ref_centered
    return float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))


def protein_axis_angle(coords):
    centered = coords - np.mean(coords, axis=0)
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    cosine = np.clip(abs(float(vt[0, 2])), 0.0, 1.0)
    return float(np.degrees(np.arccos(cosine)))


def rigid_body_axis_angles(coords):
    reference = coords[0]
    centered = reference - np.mean(reference, axis=0)
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    reference_axis = vt[0]
    angles = []
    for frame in coords:
        rotation_to_reference = kabsch_rotation(frame, reference)
        current_axis = reference_axis @ rotation_to_reference.T
        cosine = np.clip(abs(float(current_axis[2])), 0.0, 1.0)
        angles.append(np.degrees(np.arccos(cosine)))
    return np.asarray(angles)


def unwrap_xy(pos_xy, box_xy):
    box_xy = np.asarray(box_xy, dtype=np.float64)
    if box_xy.ndim == 1:
        box_xy = np.repeat(box_xy[None, :], pos_xy.shape[0], axis=0)
    if box_xy.shape != (pos_xy.shape[0], 2):
        raise ValueError("box_xy must have shape 2 or n_frame x 2")
    delta = np.diff(pos_xy, axis=0)
    step_box = box_xy[1:, None, :]
    delta -= step_box * np.round(delta / step_box)
    return np.concatenate(
        (pos_xy[:1], pos_xy[:1] + np.cumsum(delta, axis=0)), axis=0
    )


def lateral_diffusion(time, pos_xy, fit_min=1.0, fit_max=5.0):
    centered = pos_xy - np.mean(pos_xy, axis=1, keepdims=True)
    lag_time = time - time[0]
    lags = np.where((lag_time >= fit_min) & (lag_time <= fit_max))[0]
    if lags.size < 3:
        raise ValueError("trajectory is too short for the 1-5 tu diffusion fit")
    msd = np.asarray(
        [
            np.mean(np.sum((centered[lag:] - centered[:-lag]) ** 2, axis=2))
            for lag in lags
        ]
    )
    fit_time = lag_time[lags]
    slope, r_squared = linear_fit(fit_time, msd)
    alpha, _ = linear_fit(np.log(fit_time), np.log(msd))
    return slope / 4.0, r_squared, alpha


def read_segment(path):
    with h5py.File(path, "r") as h5:
        data = {
            "time": np.asarray(h5["output/time"][:], dtype=np.float64),
            "pos": np.asarray(h5["output/pos"][:, 0], dtype=np.float64),
            "hbond": np.asarray(h5["output/hbond"][:], dtype=np.float64),
            "orientation": np.asarray(
                h5["output/cgl_orientation"][:, 0], dtype=np.float64
            ),
        }
        if "output/cgl_compaction" in h5:
            data["compaction"] = np.asarray(
                h5["output/cgl_compaction"][:, 0], dtype=np.float64
            )
        if "output/mom" in h5:
            data["mom"] = np.asarray(h5["output/mom"][:, 0], dtype=np.float64)
        if "output/box" in h5:
            data["box"] = np.asarray(h5["output/box"][:], dtype=np.float64)
    return data


def concatenate_segments(paths):
    first = read_segment(paths[0])
    combined = {name: value.copy() for name, value in first.items()}
    boundary_errors = []
    for path in paths[1:]:
        segment = read_segment(path)
        errors = {
            name: float(np.max(np.abs(combined[name][-1] - segment[name][0])))
            for name in ("pos", "orientation", "compaction", "mom")
            if name in combined and name in segment
        }
        boundary_errors.append(errors)
        offset = float(combined["time"][-1])
        combined["time"] = np.concatenate(
            (combined["time"], offset + segment["time"][1:])
        )
        for name in combined:
            if name == "time":
                continue
            combined[name] = np.concatenate((combined[name], segment[name][1:]), axis=0)
    return combined, boundary_errors


def dssp_metrics(reference_pdb, positions, backbone_indices):
    structure = md.load(str(reference_pdb))
    selection = structure.topology.select(
        "protein and (name N or name CA or name C or name O)"
    )
    if selection.size != backbone_indices.size:
        raise ValueError(
            f"reference backbone has {selection.size} atoms; runtime has "
            f"{backbone_indices.size}"
        )
    topology = structure.atom_slice(selection).topology
    trajectory = md.Trajectory(positions[:, backbone_indices] / 10.0, topology)
    dssp = md.compute_dssp(trajectory, simplified=True)
    match = np.mean(dssp == dssp[0], axis=1)
    helix = np.mean(dssp == "H", axis=1)
    return match, helix, "".join(dssp[0]), "".join(dssp[-1])


def analyze(paths, reference_pdb):
    trajectory, boundary_errors = concatenate_segments(paths)
    with h5py.File(paths[0], "r") as h5:
        bb_map = h5["input/hybrid_bb_map"]
        backbone_indices = np.asarray(bb_map["atom_indices"][:], dtype=int).reshape(-1)
        ca_indices = np.asarray(bb_map["atom_indices"][:, 1], dtype=int)
        chains = decode(bb_map["bb_chain_id"][:])
        cgl_indices = np.asarray(
            h5["input/potential/compose_vector6d/elem_index"][:], dtype=int
        )
        pair = h5["input/potential/cg_lipid_pair"]
        box = np.asarray(
            [pair.attrs["x_len"], pair.attrs["y_len"], pair.attrs["z_len"]],
            dtype=np.float64,
        )
        q_group = h5.get("input/potential/cgl_response_state")
        if q_group is None:
            q_group = h5.get("input/potential/cgl_compaction_state")
        q_min = float(q_group.attrs["self_coord_min_ang"]) if q_group else None
        q_max = float(q_group.attrs["self_coord_max_ang"]) if q_group else None
        mass_scale = float(h5["input/mass"].attrs["cgl_mass_scale"])
        pair_source = str(pair.attrs.get("pair_relaxation_correction_source", ""))
        stage_dt = integration_stage_time_step(h5)

    time = trajectory["time"]
    pos = trajectory["pos"]
    ca = pos[:, ca_indices]
    rmsd = np.asarray([kabsch_rmsd(frame, ca[0]) for frame in ca])
    hbond = np.sum(trajectory["hbond"], axis=1)
    match, helix, initial_dssp, final_dssp = dssp_metrics(
        reference_pdb, pos, backbone_indices
    )

    cgl_pos = pos[:, cgl_indices]
    initial_upper = cgl_pos[0, :, 2] >= np.median(cgl_pos[0, :, 2])
    current_upper = cgl_pos[:, :, 2] >= np.median(cgl_pos[:, :, 2], axis=1)[:, None]
    flip_fraction = np.mean(current_upper != initial_upper[None, :], axis=1)
    separation = (
        np.mean(cgl_pos[:, initial_upper, 2], axis=1)
        - np.mean(cgl_pos[:, ~initial_upper, 2], axis=1)
    )
    orientation = trajectory["orientation"]
    abs_nz = np.abs(orientation[:, :, 2]) / np.maximum(
        np.linalg.norm(orientation, axis=2), 1.0e-12
    )
    box_series = trajectory.get("box")
    if box_series is None:
        box_series = np.repeat(box[None, :], time.size, axis=0)
    unwrapped_xy = unwrap_xy(cgl_pos[:, :, :2], box_series[:, :2])
    diffusion, diffusion_r2, diffusion_alpha = lateral_diffusion(time, unwrapped_xy)
    uniformity = [
        lateral_uniformity(cgl_pos, box_series, mask)
        for mask in (initial_upper, ~initial_upper)
    ]

    q = trajectory.get("compaction")
    near_values = []
    far_values = []
    if q is not None:
        protein_pos = pos[:, backbone_indices]
        for frame_cgl, frame_protein, frame_q in zip(cgl_pos, protein_pos, q):
            delta = frame_cgl[:, None, :] - frame_protein[None, :, :]
            delta -= box.reshape(1, 1, 3) * np.round(
                delta / box.reshape(1, 1, 3)
            )
            min_distance = np.sqrt(np.min(np.sum(delta * delta, axis=2), axis=1))
            near_values.extend(frame_q[min_distance <= 15.0])
            far_values.extend(frame_q[min_distance >= 25.0])

    print("segments=" + ",".join(str(path) for path in paths))
    print(
        f"  contract stage_dt={stage_dt:.6g} tu mass_scale={mass_scale:.6g} "
        f"pair_source={pair_source}"
    )
    for index, errors in enumerate(boundary_errors, start=1):
        print(
            f"  restart_boundary_{index} "
            + " ".join(f"{name}_max_error={value:.6g}" for name, value in errors.items())
        )
    print(
        f"  hbond start/min/mean/late/final={hbond[0]:.3f}/{np.min(hbond):.3f}/"
        f"{np.mean(hbond):.3f}/{np.mean(hbond[len(hbond)//2:]):.3f}/{hbond[-1]:.3f}"
    )
    print(
        f"  ca_rmsd median/late/final/max={np.median(rmsd):.3f}/"
        f"{np.mean(rmsd[len(rmsd)//2:]):.3f}/{rmsd[-1]:.3f}/{np.max(rmsd):.3f} A"
    )
    print(
        f"  dssp_match mean/late/final={np.mean(match):.3f}/"
        f"{np.mean(match[len(match)//2:]):.3f}/{match[-1]:.3f} "
        f"helix start/late/final={helix[0]:.3f}/"
        f"{np.mean(helix[len(helix)//2:]):.3f}/{helix[-1]:.3f}"
    )
    print(f"  dssp initial={initial_dssp} final={final_dssp}")
    for chain in np.unique(chains):
        chain_mask = chains == chain
        chain_ca = ca[:, chain_mask]
        angles = rigid_body_axis_angles(chain_ca)
        slope, angle_r2 = linear_fit(time, angles)
        print(
            f"  angle chain={chain} start/mean/late/final="
            f"{angles[0]:.3f}/{np.mean(angles):.3f}/"
            f"{np.mean(angles[len(angles)//2:]):.3f}/{angles[-1]:.3f} deg "
            f"range=[{np.min(angles):.3f},{np.max(angles):.3f}] "
            f"slope={slope:.6g} deg/tu R2={angle_r2:.3f}"
        )
        shape_angles = np.asarray([protein_axis_angle(frame) for frame in chain_ca])
        print(
            f"  shape_axis chain={chain} start/final/range="
            f"{shape_angles[0]:.3f}/{shape_angles[-1]:.3f}/"
            f"[{np.min(shape_angles):.3f},{np.max(shape_angles):.3f}] deg"
        )
    print(
        f"  cgl D={diffusion:.6g} A^2/tu D_x56={56.0*diffusion:.6g} "
        f"R2={diffusion_r2:.3f} alpha={diffusion_alpha:.3f} fit=[1,5] tu"
    )
    print(
        f"  leaflet start/mean/late/final={separation[0]:.3f}/"
        f"{np.mean(separation):.3f}/{np.mean(separation[len(separation)//2:]):.3f}/"
        f"{separation[-1]:.3f} A abs_nz_mean/p05={np.mean(abs_nz):.4f}/"
        f"{np.quantile(abs_nz, 0.05):.4f}"
    )
    for label, metrics in zip(("upper", "lower"), uniformity):
        print(
            f"  uniformity {label} low_k_mean/max="
            f"{np.mean(metrics['structure']):.4f}/{np.max(metrics['structure']):.4f} "
            f"occupancy6_cv_mean/max={np.mean(metrics['occupancy_cv']):.4f}/"
            f"{np.max(metrics['occupancy_cv']):.4f} "
            f"empty6_mean/max={np.mean(metrics['empty_fraction']):.4f}/"
            f"{np.max(metrics['empty_fraction']):.4f} "
            f"nearest_p05/mean={np.quantile(metrics['nearest'], 0.05):.3f}/"
            f"{np.mean(metrics['nearest']):.3f} A"
        )
    print(
        f"  leaflet_flip_fraction mean/max/final={np.mean(flip_fraction):.6f}/"
        f"{np.max(flip_fraction):.6f}/{flip_fraction[-1]:.6f}"
    )
    if q is not None:
        edge = 0.01 * (q_max - q_min)
        near_mean = float(np.mean(near_values)) if near_values else float("nan")
        far_mean = float(np.mean(far_values)) if far_values else float("nan")
        print(
            f"  compaction mean/std/near/far={np.mean(q):.4f}/{np.std(q):.4f}/"
            f"{near_mean:.4f}/{far_mean:.4f} edge_fraction="
            f"{np.mean((q <= q_min + edge) | (q >= q_max - edge)):.6f}"
        )
    else:
        print("  compaction=not_present")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("trajectory", nargs="+", type=Path)
    parser.add_argument("--reference-pdb", required=True, type=Path)
    args = parser.parse_args()
    analyze(args.trajectory, args.reference_pdb)


if __name__ == "__main__":
    main()
