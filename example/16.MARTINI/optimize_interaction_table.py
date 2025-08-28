#!/usr/bin/env python3

import sys, os
import numpy as np
import h5py
from collections import defaultdict

def optimize_interaction_table(input_file, output_file):
    """Optimize the interaction table by creating a reference-based system"""
    
    print("=" * 80)
    print("OPTIMIZING INTERACTION TABLE")
    print("=" * 80)
    
    with h5py.File(input_file, 'r') as src, h5py.File(output_file, 'w') as dst:
        # Create input group first
        dst.create_group('input')
        
        # Copy all data except the martini_potential group
        print("Copying basic data...")
        
        # Copy input data manually
        for key in src['input']:
            if key != 'potential':
                if isinstance(src['input'][key], h5py.Dataset):
                    dst['input'].create_dataset(key, data=src['input'][key][:])
                else:
                    # Copy group
                    src['input'][key].copy(src['input'][key], dst['input'])
        
        # Create potential group
        dst.create_group('input/potential')
        
        # Copy bonded interactions
        print("Copying bonded interactions...")
        for bonded_type in ['dist_spring', 'angle_spring']:
            if bonded_type in src['input/potential']:
                src['input/potential'][bonded_type].copy(src['input/potential'][bonded_type], dst['input/potential'])
        
        # Copy periodic boundary
        if 'periodic_boundary_potential' in src['input/potential']:
            src['input/potential/periodic_boundary_potential'].copy(src['input/potential/periodic_boundary_potential'], dst['input/potential'])
        
        # Optimize martini potential
        print("Optimizing martini potential...")
        martini_src = src['input/potential/martini_potential']
        martini_dst = dst.create_group('input/potential/martini_potential')
        
        # Copy attributes
        for attr_name, attr_value in martini_src.attrs.items():
            martini_dst.attrs[attr_name] = attr_value
        
        # Copy basic data
        martini_dst.create_dataset('atom_indices', data=martini_src['atom_indices'][:])
        martini_dst.create_dataset('charges', data=martini_src['charges'][:])
        
        # Optimize interaction table
        pairs = martini_src['pairs'][:]
        coefficients = martini_src['coefficients'][:]
        
        print(f"Original pairs: {len(pairs)}")
        print(f"Original coefficients: {len(coefficients)}")
        
        # Create unique coefficients table
        unique_coefficients = []
        coefficient_to_index = {}
        pair_to_coefficient_index = []
        
        for i, coeff in enumerate(coefficients):
            coeff_tuple = tuple(coeff)
            if coeff_tuple not in coefficient_to_index:
                coefficient_to_index[coeff_tuple] = len(unique_coefficients)
                unique_coefficients.append(coeff)
            pair_to_coefficient_index.append(coefficient_to_index[coeff_tuple])
        
        # Convert to numpy arrays
        unique_coefficients = np.array(unique_coefficients)
        pair_to_coefficient_index = np.array(pair_to_coefficient_index)
        
        print(f"Unique coefficients: {len(unique_coefficients)}")
        print(f"Compression ratio: {len(coefficients)} / {len(unique_coefficients)} = {len(coefficients)/len(unique_coefficients):.1f}x")
        
        # Save optimized data
        martini_dst.create_dataset('pairs', data=pairs)
        martini_dst.create_dataset('coefficients', data=unique_coefficients)
        martini_dst.create_dataset('coefficient_indices', data=pair_to_coefficient_index)
        
        # Set the optimized_format attribute to tell C++ code to use optimized format
        martini_dst.attrs['optimized_format'] = 1
        
        # Copy output data if it exists
        if 'output' in src:
            print("Copying output data...")
            src['output'].copy(src['output'], dst)
        
        # Calculate size savings
        original_size = len(coefficients) * 4 * 8  # 4 floats * 8 bytes
        optimized_size = len(unique_coefficients) * 4 * 8 + len(pair_to_coefficient_index) * 4  # unique coeffs + indices
        savings = (original_size - optimized_size) / original_size * 100
        
        print(f"\nSize optimization results:")
        print(f"  Original size: {original_size / 1024 / 1024:.1f} MB")
        print(f"  Optimized size: {optimized_size / 1024 / 1024:.1f} MB")
        print(f"  Savings: {savings:.1f}%")
        
        # Verify the optimization
        print(f"\nVerification:")
        print(f"  Pairs shape: {pairs.shape}")
        print(f"  Unique coefficients shape: {unique_coefficients.shape}")
        print(f"  Coefficient indices shape: {pair_to_coefficient_index.shape}")
        print(f"  Max coefficient index: {np.max(pair_to_coefficient_index)}")
        print(f"  Min coefficient index: {np.min(pair_to_coefficient_index)}")
        
        # Test reconstruction
        reconstructed_coefficients = unique_coefficients[pair_to_coefficient_index]
        is_correct = np.allclose(coefficients, reconstructed_coefficients)
        print(f"  Reconstruction test: {'PASSED' if is_correct else 'FAILED'}")
        
        print("\n" + "=" * 80)
        print("OPTIMIZATION COMPLETE")
        print("=" * 80)

def main():
    if len(sys.argv) != 3:
        print("Usage: python optimize_interaction_table.py <input.up> <output.up>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    if not os.path.exists(input_file):
        print(f"Error: Input file {input_file} not found")
        sys.exit(1)
    
    optimize_interaction_table(input_file, output_file)

if __name__ == "__main__":
    main()
