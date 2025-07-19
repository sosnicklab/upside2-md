#!/usr/bin/env python3

# Simple test to verify bond exclusion logic with corrected DOPC topology
bonds_list = [
    [0, 1], [0, 2], [0, 3],  # NC3 bonded to PO4, GL1, GL2
    [1, 2], [2, 3], [2, 4],  # PO4-GL1, GL1-GL2, GL1-C1A
    [4, 5], [5, 6], [6, 7],  # C1A-D2A, D2A-C3A, C3A-C4A
    [3, 8], [8, 9], [9, 10], [10, 11]  # GL2-C1B, C1B-D2B, D2B-C3B, C3B-C4B
]

# Create set of bonded pairs for 1-2 exclusions
bonded_pairs_12 = set()  # Directly bonded (1-2)

# Add 1-2 exclusions from bond list
for bond in bonds_list:
    sorted_bond = (min(bond[0], bond[1]), max(bond[0], bond[1]))
    bonded_pairs_12.add(sorted_bond)

print(f"Bonds list: {bonds_list}")
print(f"Bonded pairs set: {bonded_pairs_12}")

# Test some pairs
test_pairs = [(0, 1), (0, 2), (1, 2), (0, 3), (1, 3)]
for pair in test_pairs:
    sorted_pair = (min(pair[0], pair[1]), max(pair[0], pair[1]))
    if sorted_pair in bonded_pairs_12:
        print(f"Pair {pair} -> {sorted_pair} is EXCLUDED")
    else:
        print(f"Pair {pair} -> {sorted_pair} is INCLUDED")

# Test the specific pairs from your output
problem_pairs = [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5)]
print("\nTesting problem pairs:")
for pair in problem_pairs:
    sorted_pair = (min(pair[0], pair[1]), max(pair[0], pair[1]))
    if sorted_pair in bonded_pairs_12:
        print(f"Pair {pair} -> {sorted_pair} is EXCLUDED ✓")
    else:
        print(f"Pair {pair} -> {sorted_pair} is INCLUDED ✗") 