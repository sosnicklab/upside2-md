#!/bin/bash
set -e
set -o pipefail

RUN_DIR="$(cd "$(dirname "$0")" && pwd)"
VENV_ACTIVATE="/Users/yinhan/Documents/upside2-md/.venv/bin/activate"
if [ -f "$VENV_ACTIVATE" ]; then
    # shellcheck disable=SC1090
    source "$VENV_ACTIVATE"
fi

export TEMPERATURE="0.920"
export THERMOSTAT_TIMESCALE="0.190"
export WORK_DIR="/Users/yinhan/Documents/upside2-md/example/16.MARTINI/water_diffusion/T0p920_tau0p190"
export RUN_DIR="/Users/yinhan/Documents/upside2-md/example/16.MARTINI/water_diffusion/T0p920_tau0p190"
export CHECKPOINT_DIR="/Users/yinhan/Documents/upside2-md/example/16.MARTINI/water_diffusion/T0p920_tau0p190/checkpoints"
export GENERATE_VTF=0

bash "/Users/yinhan/Documents/upside2-md/example/16.MARTINI/run_sim_water.sh"

python3 - <<'PY'
import os
import sys

import numpy as np
import tables


def calculate_msd(xyz):
    n_frames = xyz.shape[0]
    msd = np.zeros(n_frames)

    com = np.mean(xyz, axis=1)
    xyz_centered = xyz - com[:, np.newaxis, :]

    for tau in range(1, n_frames):
        diff = xyz_centered[tau:] - xyz_centered[:-tau]
        sq_dist = np.sum(diff ** 2, axis=-1)
        msd[tau] = np.mean(sq_dist)
    return msd


pdb_id = os.environ.get("PDB_ID", "water")
checkpoint_dir = os.environ.get("CHECKPOINT_DIR", "checkpoints")
traj_file = os.path.join(checkpoint_dir, f"{pdb_id}.npt_prod.up")

if not os.path.exists(traj_file):
    print(f"ERROR: production file not found: {traj_file}")
    sys.exit(1)

try:
    with tables.open_file(traj_file, "r") as f:
        if hasattr(f.root.output, "pos"):
            xyz_raw = f.root.output.pos.read()
        elif hasattr(f.root.output, "xyz"):
            xyz_raw = f.root.output.xyz.read()
        else:
            raise ValueError("No coordinates found (pos/xyz).")

        if xyz_raw.ndim == 4:
            xyz = xyz_raw[:, 0, :, :]
        elif xyz_raw.ndim == 3:
            xyz = xyz_raw
        else:
            raise ValueError(f"Unexpected coordinate shape: {xyz_raw.shape}")

        if hasattr(f.root.output, "time"):
            time = f.root.output.time.read()
        else:
            time = np.arange(xyz.shape[0])
except Exception as exc:
    print(f"CRITICAL ERROR loading trajectory: {exc}")
    sys.exit(1)

n_frames = xyz.shape[0]
if n_frames < 4:
    print("Not enough frames for diffusion analysis.")
    sys.exit(1)

start_fraction = 0.5
start_frame = int(n_frames * start_fraction)
if start_frame >= n_frames - 2:
    print("Not enough frames after applying last-50% filter.")
    sys.exit(1)

xyz = xyz[start_frame:]
time = time[start_frame:]
if time.size > 0:
    time = time - time[0]

msd = calculate_msd(xyz)

fit_time = time
fit_msd = msd
slope, intercept = np.polyfit(fit_time, fit_msd, 1)

D_nm2ps = slope / 6.0

output_path = os.path.join(os.getcwd(), "diffusion_rate.txt")
with open(output_path, "w", encoding="utf-8") as f:
    f.write(f"Temperature: {float(os.environ.get('TEMPERATURE', 'nan')):.3f} reduced units\n")
    f.write(
        f"ThermostatTimescale: {float(os.environ.get('THERMOSTAT_TIMESCALE', 'nan')):.3f} reduced units\n"
    )
    f.write(f"DiffusionRate: {D_nm2ps:.12e} nm^2/ps\n")
    f.write(f"MSDSlope: {slope:.12e} nm^2/ps\n")
    f.write(f"AnalysisStartFrame: {start_frame}\n")
    f.write(f"AnalysisFrames: {xyz.shape[0]}\n")

print(
    f"Completed T={float(os.environ.get('TEMPERATURE', 'nan')):.3f}, "
    f"tau={float(os.environ.get('THERMOSTAT_TIMESCALE', 'nan')):.3f}: "
    f"D={D_nm2ps:.12e} nm^2/ps (Used last {xyz.shape[0]} frames)"
)
PY
