#!/usr/bin/env python3

import argparse
import os
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent
DEFAULT_BASE_DIR = SCRIPT_DIR / "lipid_diffusion"
DEFAULT_RUN_SCRIPT = SCRIPT_DIR / "run_sim_bilayer.sh"
DEFAULT_VENV_ACTIVATE = REPO_ROOT / ".venv" / "bin" / "activate"
DEFAULT_PDB_ID = "bilayer"
DEFAULT_TIME_LIMIT = "36:00:00"
DEFAULT_CPUS = 1
DEFAULT_PARTITION = ""

# Optional runtime overrides for convenience. Set to None to keep workflow defaults.
DEFAULT_EQ_NSTEPS = None
DEFAULT_PROD_NSTEPS = None
DEFAULT_EQ_TIME_STEP = None
DEFAULT_PROD_TIME_STEP = None
DEFAULT_EQ_FRAME_STEPS = None
DEFAULT_PROD_FRAME_STEPS = None


def build_temperature_range():
    return np.arange(0.600, 1.001, 0.02)


def build_tau_range():
    taus_1 = np.append(np.arange(1, 5.1, 0.5), 0.135)
    return np.sort(np.append(np.arange(0.0, 0.21, 0.01), taus_1))


def format_folder_value(value):
    return f"{value:.3f}".replace(".", "p")


def build_run_script_content(
    run_sim_bilayer,
    script_dir,
    venv_activate,
    temperature,
    tau,
    run_dir,
    pdb_id,
    eq_nsteps,
    prod_nsteps,
    eq_time_step,
    prod_time_step,
    eq_frame_steps,
    prod_frame_steps,
):
    checkpoint_dir = run_dir / "checkpoints"
    log_dir = run_dir / "logs"
    env_lines = [
        f'export TEMPERATURE="{temperature:.3f}"',
        f'export THERMOSTAT_TIMESCALE="{tau:.3f}"',
        f'export RUN_DIR="{run_dir}"',
        f'export CHECKPOINT_DIR="{checkpoint_dir}"',
        f'export LOG_DIR="{log_dir}"',
        # Keep production as NVT by default for implicit dry bilayer.
        'export PROD_70_NPT_ENABLE="${PROD_70_NPT_ENABLE:-0}"',
        f'export PDB_ID="{pdb_id}"',
    ]
    if eq_nsteps is not None:
        env_lines.extend(
            [
                f'export EQ_62_NSTEPS="{eq_nsteps}"',
                f'export EQ_63_NSTEPS="{eq_nsteps}"',
                f'export EQ_64_NSTEPS="{eq_nsteps}"',
                f'export EQ_65_NSTEPS="{eq_nsteps}"',
                f'export EQ_66_NSTEPS="{eq_nsteps}"',
            ]
        )
    if prod_nsteps is not None:
        env_lines.append(f'export PROD_70_NSTEPS="{prod_nsteps}"')
    if eq_time_step is not None:
        env_lines.append(f'export EQ_TIME_STEP="{eq_time_step}"')
    if prod_time_step is not None:
        env_lines.append(f'export PROD_TIME_STEP="{prod_time_step}"')
    if eq_frame_steps is not None:
        env_lines.append(f'export EQ_FRAME_STEPS="{eq_frame_steps}"')
    if prod_frame_steps is not None:
        env_lines.append(f'export PROD_FRAME_STEPS="{prod_frame_steps}"')

    env_block = "\n".join(env_lines)

    return f"""#!/bin/bash
set -e
set -o pipefail

RESULT_DIR="$(cd "$(dirname "$0")" && pwd)"
export RESULT_DIR
VENV_ACTIVATE="{venv_activate}"
if [ -f "$VENV_ACTIVATE" ]; then
    # shellcheck disable=SC1090
    source "$VENV_ACTIVATE"
fi

{env_block}

cd "{script_dir}"
bash "{run_sim_bilayer}" PDB_ID="$PDB_ID"

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
            raise ValueError(f"Unexpected output/pos shape: {{pos_raw.shape}}")

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

    return {{
        "slope": slope,
        "intercept": intercept,
        "d_main_a2_per_tu": d_main,
        "d_first_a2_per_tu": d_first,
        "d_second_a2_per_tu": d_second,
        "d_err_a2_per_tu": d_err,
        "fit_points": fit_t.size,
        "trim_points_each_end": trim,
    }}


def main():
    pdb_id = os.environ.get("PDB_ID", "bilayer")
    checkpoint_dir = os.environ.get("CHECKPOINT_DIR", "outputs/martini_test/checkpoints")
    result_dir = os.environ.get("RESULT_DIR", os.getcwd())
    traj_file = os.path.join(checkpoint_dir, f"{{pdb_id}}.stage_7.0.up")

    if not os.path.exists(traj_file):
        print(f"ERROR: production trajectory not found: {{traj_file}}")
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
        f.write(f"Temperature: {{float(os.environ.get('TEMPERATURE', 'nan')):.3f}} reduced_units\\n")
        f.write(f"ThermostatTimescale: {{float(os.environ.get('THERMOSTAT_TIMESCALE', 'nan')):.3f}} reduced_units\\n")
        f.write(f"PO4Beads: {{po4_idx.size}}\\n")
        f.write(f"Frames: {{pos.shape[0]}}\\n")
        f.write(f"FitPoints: {{fit['fit_points']}}\\n")
        f.write(f"TrimPointsEachEnd: {{fit['trim_points_each_end']}}\\n")
        f.write(f"MSDSlope_A2_per_tu: {{fit['slope']:.12e}}\\n")
        f.write(f"MSDIntercept_A2: {{fit['intercept']:.12e}}\\n")
        f.write(f"DiffusionRate_A2_per_tu: {{fit['d_main_a2_per_tu']:.12e}}\\n")
        f.write(f"DiffusionRate_nm2_per_tu: {{d_main_nm2_per_tu:.12e}}\\n")
        f.write(f"DiffusionHalf1_A2_per_tu: {{fit['d_first_a2_per_tu']:.12e}}\\n")
        f.write(f"DiffusionHalf2_A2_per_tu: {{fit['d_second_a2_per_tu']:.12e}}\\n")
        f.write(f"FitHalfDiff_A2_per_tu: {{fit['d_err_a2_per_tu']:.12e}}\\n")
        f.write(f"FitHalfDiff_nm2_per_tu: {{d_err_nm2_per_tu:.12e}}\\n")

    print(
        f"Completed T={{float(os.environ.get('TEMPERATURE', 'nan')):.3f}}, "
        f"tau={{float(os.environ.get('THERMOSTAT_TIMESCALE', 'nan')):.3f}}: "
        f"D={{d_main_nm2_per_tu:.12e}} nm^2/time_unit, "
        f"half-fit-diff={{d_err_nm2_per_tu:.12e}} nm^2/time_unit"
    )


if __name__ == "__main__":
    main()
PY
"""


def write_slurm_script(slurm_path, tasks_file, num_tasks, time_limit, cpus, partition):
    partition_line = f"#SBATCH --partition={partition}\n" if partition else ""
    slurm_content = (
        "#!/bin/bash\n"
        "#SBATCH --job-name=bilayer_diff_scan\n"
        "#SBATCH --output=bilayer_scan_%A_%a.out\n"
        "#SBATCH --error=bilayer_scan_%A_%a.err\n"
        "#SBATCH --ntasks=1\n"
        f"#SBATCH --cpus-per-task={cpus}\n"
        f"#SBATCH --time={time_limit}\n"
        f"#SBATCH --array=1-{num_tasks}\n"
        f"{partition_line}\n"
        "\n"
        "module load cmake\n"
        "module load openmpi\n"
        "source /home/yinhanw/project/yinhan/upside2-md/source.sh\n"
        "source /home/yinhanw/project/yinhan/upside2-md/.venv/bin/activate\n"
        "\n"
        f'TASK_LIST="{tasks_file}"\n'
        'SCRIPT_PATH=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TASK_LIST")\n'
        "\n"
        'if [ -f "$SCRIPT_PATH" ]; then\n'
        '    echo "Running task $SLURM_ARRAY_TASK_ID: $SCRIPT_PATH"\n'
        '    "$SCRIPT_PATH"\n'
        "else\n"
        '    echo "Error: Script not found at $SCRIPT_PATH"\n'
        "    exit 1\n"
        "fi\n"
    )
    slurm_path.write_text(slurm_content, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(
        description="Generate Slurm array scripts for MARTINI lipid diffusion scans."
    )
    parser.add_argument(
        "--base-dir",
        default=str(DEFAULT_BASE_DIR),
        help="Base output directory for scan runs.",
    )
    parser.add_argument(
        "--run-script",
        default=str(DEFAULT_RUN_SCRIPT),
        help="Path to run_sim_bilayer.sh.",
    )
    parser.add_argument(
        "--time-limit",
        default=DEFAULT_TIME_LIMIT,
        help="Slurm time limit (HH:MM:SS).",
    )
    parser.add_argument(
        "--cpus",
        type=int,
        default=DEFAULT_CPUS,
        help="CPUs per task.",
    )
    parser.add_argument(
        "--partition",
        default=DEFAULT_PARTITION,
        help="Optional Slurm partition name.",
    )
    parser.add_argument(
        "--pdb-id",
        default=DEFAULT_PDB_ID,
        help="PDB/system identifier passed to run_sim_bilayer.sh.",
    )
    parser.add_argument(
        "--eq-nsteps",
        type=int,
        default=DEFAULT_EQ_NSTEPS,
        help="Override equilibration stage nsteps.",
    )
    parser.add_argument(
        "--prod-nsteps",
        type=int,
        default=DEFAULT_PROD_NSTEPS,
        help="Override production stage nsteps.",
    )
    parser.add_argument(
        "--eq-time-step",
        type=float,
        default=DEFAULT_EQ_TIME_STEP,
        help="Override equilibration time step.",
    )
    parser.add_argument(
        "--prod-time-step",
        type=float,
        default=DEFAULT_PROD_TIME_STEP,
        help="Override production time step.",
    )
    parser.add_argument(
        "--eq-frame-steps",
        type=int,
        default=DEFAULT_EQ_FRAME_STEPS,
        help="Override equilibration frame-write interval.",
    )
    parser.add_argument(
        "--prod-frame-steps",
        type=int,
        default=DEFAULT_PROD_FRAME_STEPS,
        help="Override production frame-write interval.",
    )

    args = parser.parse_args()

    base_dir = Path(args.base_dir).expanduser().resolve()
    run_sim_bilayer = Path(args.run_script).expanduser().resolve()
    venv_activate = Path(os.environ.get("VENV_ACTIVATE", DEFAULT_VENV_ACTIVATE)).resolve()

    if not run_sim_bilayer.exists():
        raise FileNotFoundError(f"run_sim_bilayer.sh not found at {run_sim_bilayer}")

    temperatures = build_temperature_range()
    taus = build_tau_range()

    base_dir.mkdir(parents=True, exist_ok=True)

    run_scripts = []
    for temp in temperatures:
        for tau in taus:
            folder = f"T{format_folder_value(temp)}_tau{format_folder_value(tau)}"
            run_dir = base_dir / folder
            run_dir.mkdir(parents=True, exist_ok=True)

            script_path = run_dir / "run_simulation.sh"
            script_content = build_run_script_content(
                run_sim_bilayer=run_sim_bilayer,
                script_dir=SCRIPT_DIR,
                venv_activate=venv_activate,
                temperature=float(np.round(temp, 3)),
                tau=float(np.round(tau, 3)),
                run_dir=run_dir,
                pdb_id=args.pdb_id,
                eq_nsteps=args.eq_nsteps,
                prod_nsteps=args.prod_nsteps,
                eq_time_step=args.eq_time_step,
                prod_time_step=args.prod_time_step,
                eq_frame_steps=args.eq_frame_steps,
                prod_frame_steps=args.prod_frame_steps,
            )
            script_path.write_text(script_content, encoding="utf-8")
            script_path.chmod(0o755)
            run_scripts.append(script_path)

    tasks_file = base_dir / "tasks.txt"
    tasks_file.write_text("\n".join(str(p) for p in run_scripts) + "\n", encoding="utf-8")

    slurm_script = base_dir / "run_scan.slurm"
    write_slurm_script(
        slurm_script,
        tasks_file=tasks_file,
        num_tasks=len(run_scripts),
        time_limit=args.time_limit,
        cpus=args.cpus,
        partition=args.partition,
    )

    print(f"Generated {len(run_scripts)} simulation scripts.")
    print(f"Task list saved to: {tasks_file}")
    print(f"Master Slurm script: {slurm_script}")


if __name__ == "__main__":
    main()
