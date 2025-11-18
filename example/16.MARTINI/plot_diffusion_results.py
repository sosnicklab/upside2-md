#!/usr/bin/env python3
"""
Plotting Script for Water Diffusion Parameter Scan Results
"""

import os
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

RUN_PREFIX = "water_diffusion_scan"

def read_results():
    """Read all diffusion rate results from the scan directories"""
    results = []
    pattern = re.compile(r'T([0-9]+p[0-9]+)_tau([0-9]+p[0-9]+)')

    # Find all run directories
    if not os.path.exists(RUN_PREFIX):
        print(f"Error: Results directory '{RUN_PREFIX}' not found!")
        print("Please run the parameter scan first.")
        return []

    for run_dir in os.listdir(RUN_PREFIX):
        match = pattern.match(run_dir)
        if match:
            temp_str = match.group(1).replace("p", ".")
            tau_str = match.group(2).replace("p", ".")
            temperature = float(temp_str)
            tau = float(tau_str)

            # Read diffusion rate from the result file
            result_file = os.path.join(RUN_PREFIX, run_dir, "diffusion_rate.txt")
            if os.path.exists(result_file):
                d = None
                bulk_atoms = None
                total_atoms = None

                with open(result_file, "r") as f:
                    for line in f:
                        line = line.strip()
                        if line.startswith("DiffusionRate:"):
                            d = float(line.split(":")[1].strip().split()[0])
                        elif line.startswith("BulkWaterAtoms:"):
                            bulk_atoms = int(line.split(":")[1].strip())
                        elif line.startswith("TotalWaterAtoms:"):
                            total_atoms = int(line.split(":")[1].strip())

                if d is not None:
                    results.append((temperature, tau, d, bulk_atoms, total_atoms))

    return results

def plot_temperature_vs_diffusion(results):
    """Plot temperature vs. diffusion rate and fit a line"""
    if not results:
        return

    temperatures, taus, diff_rates, _, _ = zip(*results)

    # Create scatter plot
    plt.figure(figsize=(10, 6))
    plt.scatter(temperatures, diff_rates, c='blue', alpha=0.7, label='Data')

    # Fit linear regression
    slope, intercept, r_value, p_value, std_err = linregress(temperatures, diff_rates)
    fit_line = slope * np.array(temperatures) + intercept

    # Plot the fit
    plt.plot(temperatures, fit_line, c='red', lw=2, label=f'Fit: D = {slope:.6f}T + {intercept:.6f}\nR² = {r_value**2:.6f}')

    plt.xlabel('Temperature (reduced units)', fontsize=12)
    plt.ylabel('Diffusion Rate (cm²/s)', fontsize=12)
    plt.title('Temperature vs. Water Diffusion Rate', fontsize=14)
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    # Save plot
    plt.savefig('temperature_vs_diffusion.pdf', dpi=300, bbox_inches='tight')
    plt.savefig('temperature_vs_diffusion.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Print fit results
    print("=== TEMPERATURE vs. DIFFUSION RATE ===")
    print(f"Fit equation: D = {slope:.6f} * T + {intercept:.6f}")
    print(f"R-squared: {r_value**2:.6f}")
    print(f"Standard error: {std_err:.6f}")
    print(f"P-value: {p_value:.6f}")
    print()

def plot_thermostat_timescale_vs_diffusion(results):
    """Plot thermostat timescale vs. diffusion rate and fit a line"""
    if not results:
        return

    temperatures, taus, diff_rates, _, _ = zip(*results)

    # Create scatter plot
    plt.figure(figsize=(10, 6))
    plt.scatter(taus, diff_rates, c='green', alpha=0.7, label='Data')

    # Fit linear regression
    slope, intercept, r_value, p_value, std_err = linregress(taus, diff_rates)
    fit_line = slope * np.array(taus) + intercept

    # Plot the fit
    plt.plot(taus, fit_line, c='red', lw=2, label=f'Fit: D = {slope:.6f}τ + {intercept:.6f}\nR² = {r_value**2:.6f}')

    plt.xlabel('Thermostat Timescale (reduced units)', fontsize=12)
    plt.ylabel('Diffusion Rate (cm²/s)', fontsize=12)
    plt.title('Thermostat Timescale vs. Water Diffusion Rate', fontsize=14)
    plt.legend(fontsize=11)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    # Save plot
    plt.savefig('thermostat_timescale_vs_diffusion.pdf', dpi=300, bbox_inches='tight')
    plt.savefig('thermostat_timescale_vs_diffusion.png', dpi=300, bbox_inches='tight')
    plt.close()

    # Print fit results
    print("=== THERMOSTAT TIMESCALE vs. DIFFUSION RATE ===")
    print(f"Fit equation: D = {slope:.6f} * τ + {intercept:.6f}")
    print(f"R-squared: {r_value**2:.6f}")
    print(f"Standard error: {std_err:.6f}")
    print(f"P-value: {p_value:.6f}")
    print()

def main():
    print("=== PLOTTING DIFFUSION RATE RESULTS ===")
    print()

    # Read all results
    results = read_results()

    if not results:
        print("No results found!")
        return

    print(f"Found {len(results)} results")
    print()

    # Plot temperature vs diffusion
    plot_temperature_vs_diffusion(results)

    # Plot thermostat timescale vs diffusion
    plot_thermostat_timescale_vs_diffusion(results)

    print("=== PLOTTING COMPLETE ===")
    print("Generated plots:")
    print("  - temperature_vs_diffusion.pdf/png")
    print("  - thermostat_timescale_vs_diffusion.pdf/png")
    print()
    print("Done!")

if __name__ == "__main__":
    main()
