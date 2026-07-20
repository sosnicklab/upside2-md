#!/usr/bin/env python

import argparse
import os

import numpy as np


def load_probability(path, description):
    if not os.path.isfile(path):
        raise FileNotFoundError('cannot read {}: {}'.format(description, path))
    values = np.asarray(np.load(path), dtype=np.float64)
    if not np.all(np.isfinite(values)):
        raise ValueError('{} contains non-finite values'.format(description))
    if np.any(values < 0.0) or np.any(values > 1.0):
        raise ValueError('{} values must be in [0, 1]'.format(description))
    return values


def combine(protein_protection_path, output_path, accessibility_path=None):
    protein_protection = load_probability(protein_protection_path, 'protein protection state')
    if protein_protection.ndim != 2:
        raise ValueError('protein protection state must have shape (n_frame, n_donor)')

    if accessibility_path:
        accessibility = load_probability(accessibility_path, 'water accessibility')
        if accessibility.shape != protein_protection.shape:
            raise ValueError(
                'water accessibility shape {} does not match protein protection shape {}'.format(
                    accessibility.shape, protein_protection.shape))
        exchange_probability = (1.0 - protein_protection) * accessibility
        protection = 1.0 - exchange_probability
    else:
        protection = protein_protection

    output_directory = os.path.dirname(os.path.abspath(output_path))
    os.makedirs(output_directory, exist_ok=True)
    np.save(output_path, protection.astype(np.float32))


def main():
    parser = argparse.ArgumentParser(
        description='Combine stock Upside protein protection with an optional calibrated water-accessibility array.')
    parser.add_argument('protein_protection', help='stock protein-only PS.npy')
    parser.add_argument('output_npy', help='downstream PS.npy')
    parser.add_argument(
        '--water-accessibility',
        default=None,
        help='optional array shaped like PS with 0=inaccessible and 1=fully water-accessible',
    )
    args = parser.parse_args()
    combine(args.protein_protection, args.output_npy, args.water_accessibility)
    if args.water_accessibility:
        print('Combined protein protection and water accessibility: {}'.format(args.output_npy))
    else:
        print('No membrane accessibility correction applied: {}'.format(args.output_npy))


if __name__ == '__main__':
    main()
