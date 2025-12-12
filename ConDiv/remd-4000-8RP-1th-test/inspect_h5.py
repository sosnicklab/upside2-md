#!/usr/bin/env python3
import tables as tb
import sys
import os
import numpy as np

def analyze_file(filepath):
    """Analyzes a single H5 file for known Upside parameter tensors."""
    print(f"\n[File] {filepath}")
    
    try:
        with tb.open_file(filepath, 'r') as t:
            # 1. Check for Rotamer Interaction (sidechain.h5 or trained param file)
            if hasattr(t.root, 'pair_interaction'):
                rot_shape = t.root.pair_interaction.shape
                dim = rot_shape[-1]
                # Formula: dim = 2*n_knot_angular + 2*n_knot_sc
                # Assuming n_knot_angular is fixed at 15
                n_knot_angular = 15
                n_knot_sc = (dim - 2 * n_knot_angular) / 2
                
                print(f"  > pair_interaction shape:     {rot_shape}")
                print(f"    - Last dimension:           {dim}")
                print(f"    - Calculated n_knot_sc:     {n_knot_sc} (Default code expects 12)")

            # 2. Check for Coverage/HBond (sidechain.h5 or trained param file)
            if hasattr(t.root, 'coverage_interaction'):
                cov_shape = t.root.coverage_interaction.shape
                dim = cov_shape[-1]
                # Formula: dim = 2*n_knot_angular + 2*n_knot_hb
                n_knot_angular = 15
                n_knot_hb = (dim - 2 * n_knot_angular) / 2
                
                print(f"  > coverage_interaction shape: {cov_shape}")
                print(f"    - Last dimension:           {dim}")
                print(f"    - Calculated n_knot_hb:     {n_knot_hb} (Default code expects 10)")

            # 3. Check for Environment Energies (environment.h5)
            if hasattr(t.root, 'energies'):
                print(f"  > energies shape:             {t.root.energies.shape}")

    except Exception as e:
        print(f"  [!] Error reading file: {e}")

def main():
    if len(sys.argv) < 2:
        print("Usage: python inspect_h5.py <file.h5_or_directory>")
        sys.exit(1)
    
    target_path = sys.argv[1]
    
    if not os.path.exists(target_path):
        print(f"Error: Path '{target_path}' not found.")
        sys.exit(1)

    print(f"--- Inspecting Parameters ---")

    if os.path.isfile(target_path):
        # Case 1: User provided a specific H5 file
        analyze_file(target_path)
    
    elif os.path.isdir(target_path):
        # Case 2: User provided a directory (check for standard init files)
        sc_file = os.path.join(target_path, 'sidechain.h5')
        env_file = os.path.join(target_path, 'environment.h5')
        
        if os.path.exists(sc_file):
            analyze_file(sc_file)
        else:
            print(f"  [!] No 'sidechain.h5' found in directory.")

        if os.path.exists(env_file):
            analyze_file(env_file)
    
    print("\n------------------------------------------------")

if __name__ == "__main__":
    main()
