#!/bin/bash
set -e
set -o pipefail

RESULT_DIR="$(cd "$(dirname "$0")" && pwd)"
export RESULT_DIR
VENV_ACTIVATE="/Users/yinhan/Documents/upside2-md/.venv/bin/activate"
if [ -f "$VENV_ACTIVATE" ]; then
    # shellcheck disable=SC1090
    source "$VENV_ACTIVATE"
fi

export TEMPERATURE="0.900"
export THERMOSTAT_TIMESCALE="8.000"
export RUN_DIR="/Users/yinhan/Documents/upside2-md/example/16.MARTINI/lipid_diffusion/T0p900_tau8p000"
export CHECKPOINT_DIR="/Users/yinhan/Documents/upside2-md/example/16.MARTINI/lipid_diffusion/T0p900_tau8p000/checkpoints"
export LOG_DIR="/Users/yinhan/Documents/upside2-md/example/16.MARTINI/lipid_diffusion/T0p900_tau8p000/logs"
export PROD_70_NPT_ENABLE="${PROD_70_NPT_ENABLE:-0}"
export PDB_ID="bilayer"

cd "/Users/yinhan/Documents/upside2-md/example/16.MARTINI"
bash "/Users/yinhan/Documents/upside2-md/example/16.MARTINI/run_sim_bilayer.sh" PDB_ID="$PDB_ID"

python3 - <<'PY'
import os
import sys

import numpy as np
import tables


def _decode_bytes_array(arr):
    out = []
    for x in arr:
        if isinstance(x, (bytes, bytearray)):
            out.append(x.decode("utf-8").strip())
        else:
            out.append(str(x).strip())
    return np.array(out, dtype=object)


def _load_traj(path):
    with tables.open_file(path, "r") as h5:
        if not hasattr(h5.root.output, "pos"):
            raise ValueError("Missing output/pos trajectory dataset.")
        pos_raw = h5.root.output.pos.read()
        if pos_raw.ndim == 4:
            pos = pos_raw[:, 0, :, :]
        elif pos_raw.ndim == 3:
            pos = pos_raw
        else:
            raise ValueError(f"Unexpected output/pos shape: {pos_raw.shape}")

        if hasattr(h5.root.output, "time"):
            time = np.asarray(h5.root.output.time.read(), dtype=float)
        else:
            time = np.arange(pos.shape[0], dtype=float)

        if hasattr(h5.root.output, "box"):
            box_raw = np.asarray(h5.root.output.box.read(), dtype=float)
            if box_raw.ndim == 3:
                box = box_raw[:, 0, :]
            elif box_raw.ndim == 2:
                box = box_raw
            else:
                box = None
        else:
            box = None

        atom_names = _decode_bytes_array(h5.root.input.atom_names.read())
        masses = np.asarray(h5.root.input.mass.read(), dtype=float) if hasattr(h5.root.input, "mass") else None
    return pos.astype(float), time, box, atom_names, masses


def _unwrap_xy(pos, box):
    out = pos.copy()
    n_frames = pos.shape[0]
    if n_frames < 2:
        return out

    for i in range(1, n_frames):
        delta_xy = pos[i, :, :2] - pos[i - 1, :, :2]
        if box is not None and i - 1 < box.shape[0]:
            lx = box[i - 1, 0]
            ly = box[i - 1, 1]
            if lx > 0:
                delta_xy[:, 0] -= np.round(delta_xy[:, 0] / lx) * lx
            if ly > 0:
                delta_xy[:, 1] -= np.round(delta_xy[:, 1] / ly) * ly
        out[i, :, :2] = out[i - 1, :, :2] + delta_xy
        out[i, :, 2] = pos[i, :, 2]
    return out


def _remove_com_xy(pos, masses):
    if masses is not None and masses.size == pos.shape[1] and np.all(masses > 0):
        com_xy = np.average(pos[:, :, :2], axis=1, weights=masses)
    else:
        com_xy = np.mean(pos[:, :, :2], axis=1)
    centered = pos.copy()
    centered[:, :, :2] -= com_xy[:, np.newaxis, :]
    return centered


def _compute_lateral_msd(po4_xy, time):
    n_frames = po4_xy.shape[0]
    if n_frames < 4:
        raise ValueError("Not enough frames for MSD analysis.")

    lag_times = np.zeros(n_frames - 1, dtype=float)
    msd = np.zeros(n_frames - 1, dtype=float)

    for lag in range(1, n_frames):
        disp = po4_xy[lag:] - po4_xy[:-lag]
        sq = np.sum(disp * disp, axis=2)  # (n_pairs, n_po4)
        msd[lag - 1] = np.mean(sq)
        lag_times[lag - 1] = np.mean(time[lag:] - time[:-lag])
    return lag_times, msd


def _fit_diffusion_2d(lag_times, msd):
    n = lag_times.size
    trim = int(0.1 * n)
    lo = trim
    hi = n - trim
    if hi - lo < 4:
        raise ValueError("Insufficient points after 10% trimming for fit.")

    fit_t = lag_times[lo:hi]
    fit_msd = msd[lo:hi]
    slope, intercept = np.polyfit(fit_t, fit_msd, 1)
    d_main = slope / 4.0

    mid = fit_t.size // 2
    if mid < 2 or fit_t.size - mid < 2:
        d_first = np.nan
        d_second = np.nan
        d_err = np.nan
    else:
        slope_first, _ = np.polyfit(fit_t[:mid], fit_msd[:mid], 1)
        slope_second, _ = np.polyfit(fit_t[mid:], fit_msd[mid:], 1)
        d_first = slope_first / 4.0
        d_second = slope_second / 4.0
        d_err = abs(d_first - d_second)

    return {
        "slope": slope,
        "intercept": intercept,
        "d_main_a2_per_tu": d_main,
        "d_first_a2_per_tu": d_first,
        "d_second_a2_per_tu": d_second,
        "d_err_a2_per_tu": d_err,
        "fit_points": fit_t.size,
        "trim_points_each_end": trim,
    }


def main():
    pdb_id = os.environ.get("PDB_ID", "bilayer")
    checkpoint_dir = os.environ.get("CHECKPOINT_DIR", "outputs/martini_test/checkpoints")
    result_dir = os.environ.get("RESULT_DIR", os.getcwd())
    traj_file = os.path.join(checkpoint_dir, f"{pdb_id}.stage_7.0.up")

    if not os.path.exists(traj_file):
        print(f"ERROR: production trajectory not found: {traj_file}")
        sys.exit(1)

    pos, time, box, atom_names, masses = _load_traj(traj_file)
    po4_idx = np.where(atom_names == "PO4")[0]
    if po4_idx.size == 0:
        print("ERROR: No PO4 beads found in input/atom_names.")
        sys.exit(1)

    pos_unwrapped = _unwrap_xy(pos, box)
    pos_centered = _remove_com_xy(pos_unwrapped, masses)
    po4_xy = pos_centered[:, po4_idx, :2]

    lag_times, msd = _compute_lateral_msd(po4_xy, time)
    fit = _fit_diffusion_2d(lag_times, msd)

    # 1 A^2 = 0.01 nm^2
    d_main_nm2_per_tu = fit["d_main_a2_per_tu"] * 0.01
    d_err_nm2_per_tu = fit["d_err_a2_per_tu"] * 0.01 if np.isfinite(fit["d_err_a2_per_tu"]) else np.nan

    out_file = os.path.join(result_dir, "diffusion_rate.txt")
    with open(out_file, "w", encoding="utf-8") as f:
        f.write(f"Temperature: {float(os.environ.get('TEMPERATURE', 'nan')):.3f} reduced_units\n")
        f.write(f"ThermostatTimescale: {float(os.environ.get('THERMOSTAT_TIMESCALE', 'nan')):.3f} reduced_units\n")
        f.write(f"PO4Beads: {po4_idx.size}\n")
        f.write(f"Frames: {pos.shape[0]}\n")
        f.write(f"FitPoints: {fit['fit_points']}\n")
        f.write(f"TrimPointsEachEnd: {fit['trim_points_each_end']}\n")
        f.write(f"MSDSlope_A2_per_tu: {fit['slope']:.12e}\n")
        f.write(f"MSDIntercept_A2: {fit['intercept']:.12e}\n")
        f.write(f"DiffusionRate_A2_per_tu: {fit['d_main_a2_per_tu']:.12e}\n")
        f.write(f"DiffusionRate_nm2_per_tu: {d_main_nm2_per_tu:.12e}\n")
        f.write(f"DiffusionHalf1_A2_per_tu: {fit['d_first_a2_per_tu']:.12e}\n")
        f.write(f"DiffusionHalf2_A2_per_tu: {fit['d_second_a2_per_tu']:.12e}\n")
        f.write(f"FitHalfDiff_A2_per_tu: {fit['d_err_a2_per_tu']:.12e}\n")
        f.write(f"FitHalfDiff_nm2_per_tu: {d_err_nm2_per_tu:.12e}\n")

    print(
        f"Completed T={float(os.environ.get('TEMPERATURE', 'nan')):.3f}, "
        f"tau={float(os.environ.get('THERMOSTAT_TIMESCALE', 'nan')):.3f}: "
        f"D={d_main_nm2_per_tu:.12e} nm^2/time_unit, "
        f"half-fit-diff={d_err_nm2_per_tu:.12e} nm^2/time_unit"
    )


if __name__ == "__main__":
    main()
PY
