import argparse
import os
import re
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_BASE_DIR = SCRIPT_DIR / "water_diffusion"


def read_results(base_path):
    results = []
    pattern = re.compile(r"T([0-9]+p[0-9]+)_tau([0-9]+p[0-9]+)")

    if not os.path.exists(base_path):
        print(f"Error: Results directory '{base_path}' not found!")
        return results

    print(f"Scanning directory: {base_path}")

    for run_dir in sorted(os.listdir(base_path)):
        match = pattern.match(run_dir)
        if not match:
            continue

        temp_str = match.group(1).replace("p", ".")
        tau_str = match.group(2).replace("p", ".")
        temperature = float(temp_str)
        tau = float(tau_str)

        result_file = os.path.join(base_path, run_dir, "diffusion_rate.txt")
        if not os.path.exists(result_file):
            continue

        diffusion = None
        with open(result_file, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line.startswith("DiffusionRate:"):
                    parts = line.split(":")
                    if len(parts) > 1:
                        value_part = parts[1].strip().split()[0]
                        try:
                            diffusion = float(value_part)
                        except ValueError:
                            diffusion = None
                    break

        if diffusion is not None:
            results.append({"T": temperature, "tau": tau, "D": diffusion})
            print(f"  [Loaded] T={temperature:.3f}, tau={tau:.3f} -> D={diffusion:.3e}")

    return results


def linear_fit(x, y):
    if len(x) < 2:
        return None, None, None

    x_arr = np.array(x)
    y_arr = np.array(y)
    slope, intercept = np.polyfit(x_arr, y_arr, 1)
    y_pred = slope * x_arr + intercept

    ss_res = np.sum((y_arr - y_pred) ** 2)
    ss_tot = np.sum((y_arr - np.mean(y_arr)) ** 2)
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return slope, intercept, r2


def safe_png_name(stem):
    return stem.replace(".", "p") + ".png"


def plot_diffusion_vs_temp_fixed_tau(results, target_tau, output_dir):
    subset = [r for r in results if abs(r["tau"] - target_tau) < 1e-4]
    subset.sort(key=lambda x: x["T"])

    if not subset:
        print(f"\n[Warning] No results found for tau = {target_tau}")
        available_taus = sorted({r["tau"] for r in results})
        print(f"Available thermostat timescales found: {available_taus}")
        return

    temps = [r["T"] for r in subset]
    diffs = [r["D"] for r in subset]

    plt.figure(figsize=(8, 6))
    plt.plot(temps, diffs, "o-", linewidth=2, markersize=8, label=f"tau={target_tau}")

    slope, intercept, r2 = linear_fit(temps, diffs)
    if slope is not None:
        fit_line = [slope * t + intercept for t in temps]
        plt.plot(temps, fit_line, "--", color="red", label=f"Fit: R^2={r2:.3f}")

    plt.xlabel("Temperature (reduced units)", fontsize=12)
    plt.ylabel("Diffusion Rate (nm^2/ps)", fontsize=12)
    plt.title(f"Water Diffusion vs Temperature (tau={target_tau})", fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    output_filename = output_dir / safe_png_name(f"diffusion_vs_T_tau{target_tau:.3f}")
    plt.savefig(output_filename, dpi=300)
    print(f"\n[Success] Plot saved to {output_filename}")


def plot_diffusion_vs_tau_fixed_temp(results, target_temp, output_dir):
    subset = [r for r in results if abs(r["T"] - target_temp) < 1e-4]
    subset.sort(key=lambda x: x["tau"])

    if not subset:
        print(f"\n[Warning] No results found for T = {target_temp}")
        available_temps = sorted({r["T"] for r in results})
        print(f"Available temperatures found: {available_temps}")
        return

    taus = [r["tau"] for r in subset]
    diffs = [r["D"] for r in subset]

    plt.figure(figsize=(8, 6))
    plt.plot(taus, diffs, "o-", linewidth=2, markersize=8, label=f"T={target_temp}")

    plt.xlabel("Thermostat Timescale (reduced units)", fontsize=12)
    plt.ylabel("Diffusion Rate (nm^2/ps)", fontsize=12)
    plt.title(f"Water Diffusion vs Thermostat Timescale (T={target_temp})", fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()

    output_filename = output_dir / safe_png_name(f"diffusion_vs_tau_T{target_temp:.3f}")
    plt.savefig(output_filename, dpi=300)
    print(f"\n[Success] Plot saved to {output_filename}")


def select_target(values, provided, label):
    if provided is not None:
        return provided
    values_sorted = sorted(values)
    if not values_sorted:
        raise ValueError(f"No {label} values available to select from.")
    return values_sorted[len(values_sorted) // 2]


def main():
    parser = argparse.ArgumentParser(description="Plot MARTINI water diffusion scan results.")
    parser.add_argument(
        "--base-dir",
        default=str(DEFAULT_BASE_DIR),
        help="Base directory containing scan results.",
    )
    parser.add_argument(
        "--output-dir",
        default=None,
        help="Directory to save plots (defaults to base directory).",
    )
    parser.add_argument(
        "--target-tau",
        type=float,
        default=None,
        help="Thermostat timescale to plot diffusion vs temperature.",
    )
    parser.add_argument(
        "--target-temp",
        type=float,
        default=None,
        help="Temperature to plot diffusion vs thermostat timescale.",
    )

    args = parser.parse_args()

    base_dir = Path(args.base_dir).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve() if args.output_dir else base_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== ANALYSIS STARTED ===")
    results = read_results(str(base_dir))

    if not results:
        print("No results parsed. Check if simulations completed.")
        return

    all_taus = {r["tau"] for r in results}
    all_temps = {r["T"] for r in results}

    target_tau = select_target(all_taus, args.target_tau, "tau")
    target_temp = select_target(all_temps, args.target_temp, "temperature")

    print(f"\n--- Plotting for Target Tau = {target_tau:.3f} ---")
    plot_diffusion_vs_temp_fixed_tau(results, target_tau=target_tau, output_dir=output_dir)

    print(f"\n--- Plotting for Target Temperature = {target_temp:.3f} ---")
    plot_diffusion_vs_tau_fixed_temp(results, target_temp=target_temp, output_dir=output_dir)

    print("\n=== ANALYSIS COMPLETED ===")


if __name__ == "__main__":
    main()
