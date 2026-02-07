import argparse
import os
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent
DEFAULT_BASE_DIR = SCRIPT_DIR / "water_diffusion"
DEFAULT_RUN_SCRIPT = SCRIPT_DIR / "run_sim_water.sh"
DEFAULT_VENV_ACTIVATE = REPO_ROOT / ".venv" / "bin" / "activate"


def build_temperature_range():
    return np.arange(0.600, 1.001, 0.02)


def build_tau_range():
    taus_1 = np.append(np.arange(1, 5.1, 0.5), 0.135)
    return np.sort(np.append(np.arange(0.0, 0.21, 0.01), taus_1))


def format_folder_value(value):
    return f"{value:.3f}".replace(".", "p")


def build_run_script_content(
    run_sim_water,
    venv_activate,
    temperature,
    tau,
    run_dir,
    generate_vtf,
    npt_equil_steps,
    npt_prod_steps,
    frame_interval,
):
    vtf_flag = 1 if generate_vtf else 0
    env_lines = [
        f'export TEMPERATURE="{temperature:.3f}"',
        f'export THERMOSTAT_TIMESCALE="{tau:.3f}"',
        f'export WORK_DIR="{run_dir}"',
        f'export RUN_DIR="{run_dir}"',
        f'export CHECKPOINT_DIR="{run_dir}/checkpoints"',
        f"export GENERATE_VTF={vtf_flag}",
    ]
    if npt_equil_steps is not None:
        env_lines.append(f'export NPT_EQUIL_STEPS="{npt_equil_steps}"')
    if npt_prod_steps is not None:
        env_lines.append(f'export NPT_PROD_STEPS="{npt_prod_steps}"')
    if frame_interval is not None:
        env_lines.append(f'export FRAME_INTERVAL="{frame_interval}"')

    env_block = "\n".join(env_lines)

    return f"""#!/bin/bash
set -e

RUN_DIR="$(cd "$(dirname "$0")" && pwd)"
VENV_ACTIVATE="{venv_activate}"
if [ -f "$VENV_ACTIVATE" ]; then
    # shellcheck disable=SC1090
    source "$VENV_ACTIVATE"
fi

{env_block}

bash "{run_sim_water}"

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
traj_file = os.path.join(checkpoint_dir, f"{{pdb_id}}.npt_prod.up")

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
    f.write(f"Temperature: {float(os.environ.get('TEMPERATURE', 'nan')):.3f} reduced units\\n")
    f.write(
        f"ThermostatTimescale: {float(os.environ.get('THERMOSTAT_TIMESCALE', 'nan')):.3f} reduced units\\n"
    )
    f.write(f"DiffusionRate: {D_nm2ps:.12e} nm^2/ps\\n")
    f.write(f"MSDSlope: {slope:.12e} nm^2/ps\\n")
    f.write(f"AnalysisStartFrame: {start_frame}\\n")
    f.write(f"AnalysisFrames: {xyz.shape[0]}\\n")

print(
    f"Completed T={float(os.environ.get('TEMPERATURE', 'nan')):.3f}, "
    f"tau={float(os.environ.get('THERMOSTAT_TIMESCALE', 'nan')):.3f}: "
    f"D={D_nm2ps:.12e} nm^2/ps (Used last {xyz.shape[0]} frames)"
)
PY
"""


def write_slurm_script(slurm_path, tasks_file, num_tasks, time_limit, cpus, partition):
    partition_line = f"#SBATCH --partition={partition}\n" if partition else ""
    slurm_content = (
        "#!/bin/bash\n"
        "#SBATCH --job-name=water_diff_scan\n"
        "#SBATCH --output=water_scan_%A_%a.out\n"
        "#SBATCH --error=water_scan_%A_%a.err\n"
        "#SBATCH --ntasks=1\n"
        f"#SBATCH --cpus-per-task={cpus}\n"
        f"#SBATCH --time={time_limit}\n"
        f"#SBATCH --array=1-{num_tasks}\n"
        f"{partition_line}\n"
        f'TASK_LIST="{tasks_file}"\n'
        'SCRIPT_PATH=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TASK_LIST")\n\n'
        'if [ -f "$SCRIPT_PATH" ]; then\n'
        '    echo "Running task $SLURM_ARRAY_TASK_ID: $SCRIPT_PATH"\n'
        '    "$SCRIPT_PATH"\n'
        "else\n"
        '    echo "Error: Script not found at $SCRIPT_PATH"\n'
        "    exit 1\n"
        "fi\n"
    )

    with open(slurm_path, "w", encoding="utf-8") as f:
        f.write(slurm_content)


def main():
    parser = argparse.ArgumentParser(
        description="Generate Slurm array scripts for MARTINI water diffusion scans."
    )
    parser.add_argument(
        "--base-dir",
        default=str(DEFAULT_BASE_DIR),
        help="Base output directory for scan runs.",
    )
    parser.add_argument(
        "--run-script",
        default=str(DEFAULT_RUN_SCRIPT),
        help="Path to run_sim_water.sh.",
    )
    parser.add_argument(
        "--time-limit",
        default="36:00:00",
        help="Slurm time limit (HH:MM:SS).",
    )
    parser.add_argument(
        "--cpus",
        type=int,
        default=1,
        help="CPUs per task.",
    )
    parser.add_argument(
        "--partition",
        default="",
        help="Optional Slurm partition name.",
    )
    parser.add_argument(
        "--generate-vtf",
        action="store_true",
        help="Enable VTF generation in each run.",
    )
    parser.add_argument(
        "--npt-equil-steps",
        type=int,
        default=None,
        help="Override NPT equilibration steps.",
    )
    parser.add_argument(
        "--npt-prod-steps",
        type=int,
        default=None,
        help="Override NPT production steps.",
    )
    parser.add_argument(
        "--frame-interval",
        type=int,
        default=None,
        help="Override frame interval.",
    )

    args = parser.parse_args()

    base_dir = Path(args.base_dir).expanduser().resolve()
    run_sim_water = Path(args.run_script).expanduser().resolve()
    venv_activate = Path(os.environ.get("VENV_ACTIVATE", DEFAULT_VENV_ACTIVATE)).resolve()

    if not run_sim_water.exists():
        raise FileNotFoundError(f"run_sim_water.sh not found at {run_sim_water}")

    base_dir.mkdir(parents=True, exist_ok=True)

    temperatures = build_temperature_range()
    taus = build_tau_range()

    run_scripts = []
    for temp in temperatures:
        for tau in taus:
            folder_name = f"T{format_folder_value(temp)}_tau{format_folder_value(tau)}"
            run_dir = base_dir / folder_name
            run_dir.mkdir(parents=True, exist_ok=True)

            script_content = build_run_script_content(
                run_sim_water=run_sim_water,
                venv_activate=venv_activate,
                temperature=float(np.round(temp, 3)),
                tau=float(np.round(tau, 3)),
                run_dir=run_dir,
                generate_vtf=args.generate_vtf,
                npt_equil_steps=args.npt_equil_steps,
                npt_prod_steps=args.npt_prod_steps,
                frame_interval=args.frame_interval,
            )

            script_path = run_dir / "run_simulation.sh"
            script_path.write_text(script_content, encoding="utf-8")
            script_path.chmod(0o755)
            run_scripts.append(script_path)

    tasks_file = base_dir / "tasks.txt"
    tasks_file.write_text(
        "\n".join(str(path) for path in run_scripts) + "\n", encoding="utf-8"
    )

    slurm_script_path = base_dir / "run_scan.slurm"
    write_slurm_script(
        slurm_script_path,
        tasks_file,
        num_tasks=len(run_scripts),
        time_limit=args.time_limit,
        cpus=args.cpus,
        partition=args.partition,
    )

    print(f"Generated {len(run_scripts)} simulation scripts.")
    print(f"Task list saved to: {tasks_file}")
    print(f"Generated master SLURM script: {slurm_script_path}")


if __name__ == "__main__":
    main()
