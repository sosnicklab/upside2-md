import os
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

# --- CONFIGURATION ---
# Correct path to your simulation results
RUN_PREFIX = "/home/yinhanw/project/yinhan/upside2-md/example/16.MARTINI/water_diffusion"
TARGET_TAU = 0.135  # The specific thermostat interval you want to analyze

def read_results(base_path):
    """Read all diffusion rate results from the scan directories"""
    results = []
    # Regex to match folder names like T0p600_tau0p000
    pattern = re.compile(r'T([0-9]+p[0-9]+)_tau([0-9]+p[0-9]+)')

    if not os.path.exists(base_path):
        print(f"Error: Results directory '{base_path}' not found!")
        return []

    print(f"Scanning directory: {base_path}")
    
    # Iterate over all directories
    for run_dir in sorted(os.listdir(base_path)):
        match = pattern.match(run_dir)
        if match:
            # Parse parameters from folder name
            temp_str = match.group(1).replace("p", ".")
            tau_str = match.group(2).replace("p", ".")
            temperature = float(temp_str)
            tau = float(tau_str)

            result_file = os.path.join(base_path, run_dir, "diffusion_rate.txt")
            
            if os.path.exists(result_file):
                d = None
                with open(result_file, "r") as f:
                    for line in f:
                        line = line.strip()
                        # Robustly parse "DiffusionRate: X.XXXX cm^2/s"
                        if line.startswith("DiffusionRate:"):
                            parts = line.split(":")
                            if len(parts) > 1:
                                value_part = parts[1].strip().split()[0] # Get number before units
                                try:
                                    d = float(value_part)
                                except ValueError:
                                    pass
                
                if d is not None:
                    results.append({'T': temperature, 'tau': tau, 'D': d})
                    print(f"  [Loaded] T={temperature:.3f}, tau={tau:.3f} -> D={d:.3e}")
            else:
                 # Directory exists but job might not have finished/written result yet
                 pass

    return results

def plot_diffusion_vs_temp_fixed_tau(results, target_tau):
    """Plot Diffusion vs Temp for a specific Tau"""
    
    # Filter results for the target tau (using epsilon for float comparison)
    subset = [r for r in results if abs(r['tau'] - target_tau) < 1e-4]
    
    # Sort by temperature for clean plotting
    subset.sort(key=lambda x: x['T'])
    
    if not subset:
        print(f"\n[Warning] No results found for tau = {target_tau}")
        
        # Show what IS available to help debug
        available_taus = sorted(list(set(r['tau'] for r in results)))
        print(f"Available thermostat timescales found: {available_taus}")
        return

    temps = [r['T'] for r in subset]
    diffs = [r['D'] for r in subset]

    # Plotting
    plt.figure(figsize=(8, 6))
    plt.plot(temps, diffs, 'o-', linewidth=2, markersize=8, label=f'tau={target_tau}')
    
    # Optional: Linear Fit
    if len(temps) > 1:
        slope, intercept, r_value, _, _ = linregress(temps, diffs)
        fit_line = [slope * t + intercept for t in temps]
        plt.plot(temps, fit_line, '--', color='red', label=f'Fit: R²={r_value**2:.3f}')

    plt.xlabel('Temperature (reduced units)', fontsize=12)
    plt.ylabel('Diffusion Rate (cm²/s)', fontsize=12)
    plt.title(f'Water Diffusion vs Temperature (tau={target_tau})', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    
    # Save output
    output_filename = f'diffusion_vs_T_tau{target_tau:.3f}.png'.replace('.', 'p')
    plt.savefig(output_filename, dpi=300)
    print(f"\n[Success] Plot saved to {output_filename}")

def main():
    print("=== ANALYSIS STARTED ===")
    results = read_results(RUN_PREFIX)
    
    if not results:
        print("No results parsed. Check if simulations completed.")
        return

    # Plot specific request
    print(f"\n--- Plotting for Target Tau = {TARGET_TAU} ---")
    plot_diffusion_vs_temp_fixed_tau(results, target_tau=TARGET_TAU)

if __name__ == "__main__":
    main()