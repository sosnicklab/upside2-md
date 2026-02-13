#!/usr/bin/env python3

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_BASE_DIR = SCRIPT_DIR / "lipid_diffusion"


def read_scan_results(base_dir):
    pattern = re.compile(r"T([0-9]+p[0-9]+)_tau([0-9]+p[0-9]+)")
    results = []

    if not base_dir.exists():
        print(f"Results directory not found: {base_dir}")
        return results

    for run_dir in sorted(base_dir.iterdir()):
        if not run_dir.is_dir():
            continue
        match = pattern.match(run_dir.name)
        if not match:
            continue

        temperature = float(match.group(1).replace("p", "."))
        tau = float(match.group(2).replace("p", "."))
        result_file = run_dir / "diffusion_rate.txt"
        if not result_file.exists():
            continue

        d_nm2 = None
        d_err_nm2 = None
        with result_file.open("r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line.startswith("DiffusionRate_nm2_per_tu:"):
                    try:
                        d_nm2 = float(line.split(":", 1)[1].strip())
                    except ValueError:
                        d_nm2 = None
                elif line.startswith("FitHalfDiff_nm2_per_tu:"):
                    try:
                        d_err_nm2 = float(line.split(":", 1)[1].strip())
                    except ValueError:
                        d_err_nm2 = None

        if d_nm2 is None:
            continue
        if d_err_nm2 is None:
            d_err_nm2 = np.nan

        results.append({"T": temperature, "tau": tau, "D": d_nm2, "D_err": d_err_nm2})

    return results


def plot_diffusion_vs_temperature(results, output_dir):
    taus = sorted({r["tau"] for r in results})
    plt.figure(figsize=(9, 6))

    for tau in taus:
        subset = [r for r in results if abs(r["tau"] - tau) < 1e-6]
        subset.sort(key=lambda x: x["T"])
        if not subset:
            continue
        x = [r["T"] for r in subset]
        y = [r["D"] for r in subset]
        yerr = [r["D_err"] if np.isfinite(r["D_err"]) else 0.0 for r in subset]
        plt.errorbar(x, y, yerr=yerr, marker="o", linewidth=1.5, capsize=3, label=f"tau={tau:g}")

    plt.xlabel("Temperature (reduced units)")
    plt.ylabel("Lateral Diffusion D (nm^2 / time_unit)")
    plt.title("DOPC Lateral Diffusion vs Temperature")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    out = output_dir / "lipid_diffusion_vs_temperature.png"
    plt.savefig(out, dpi=300)
    plt.close()
    print(f"Saved: {out}")


def plot_diffusion_vs_tau(results, output_dir):
    temperatures = sorted({r["T"] for r in results})
    plt.figure(figsize=(9, 6))

    for temp in temperatures:
        subset = [r for r in results if abs(r["T"] - temp) < 1e-6]
        subset.sort(key=lambda x: x["tau"])
        if not subset:
            continue
        x = [r["tau"] for r in subset]
        y = [r["D"] for r in subset]
        yerr = [r["D_err"] if np.isfinite(r["D_err"]) else 0.0 for r in subset]
        plt.errorbar(x, y, yerr=yerr, marker="o", linewidth=1.5, capsize=3, label=f"T={temp:g}")

    plt.xlabel("Thermostat Timescale (reduced units)")
    plt.ylabel("Lateral Diffusion D (nm^2 / time_unit)")
    plt.title("DOPC Lateral Diffusion vs Thermostat Timescale")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    out = output_dir / "lipid_diffusion_vs_tau.png"
    plt.savefig(out, dpi=300)
    plt.close()
    print(f"Saved: {out}")


def plot_heatmap(results, output_dir):
    temperatures = sorted({r["T"] for r in results})
    taus = sorted({r["tau"] for r in results})
    if not temperatures or not taus:
        return

    grid = np.full((len(temperatures), len(taus)), np.nan, dtype=float)
    t_index = {t: i for i, t in enumerate(temperatures)}
    tau_index = {tau: j for j, tau in enumerate(taus)}

    for r in results:
        grid[t_index[r["T"]], tau_index[r["tau"]]] = r["D"]

    plt.figure(figsize=(8, 6))
    im = plt.imshow(
        grid,
        origin="lower",
        aspect="auto",
        cmap="viridis",
        extent=[min(taus), max(taus), min(temperatures), max(temperatures)],
    )
    plt.colorbar(im, label="D (nm^2 / time_unit)")
    plt.xlabel("Thermostat Timescale (reduced units)")
    plt.ylabel("Temperature (reduced units)")
    plt.title("DOPC Lateral Diffusion Heatmap")
    plt.tight_layout()
    out = output_dir / "lipid_diffusion_heatmap.png"
    plt.savefig(out, dpi=300)
    plt.close()
    print(f"Saved: {out}")


def main():
    base_dir = DEFAULT_BASE_DIR.resolve()
    output_dir = base_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Scanning results in: {base_dir}")
    results = read_scan_results(base_dir)
    if not results:
        print("No diffusion results found.")
        return

    print(f"Loaded {len(results)} data points.")
    plot_diffusion_vs_temperature(results, output_dir)
    plot_diffusion_vs_tau(results, output_dir)
    plot_heatmap(results, output_dir)
    print("Plot generation complete.")


if __name__ == "__main__":
    main()
