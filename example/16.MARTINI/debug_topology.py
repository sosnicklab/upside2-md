#!/usr/bin/env python3

# Simple debug script to see the .itp file structure
itp_file = "martini_v2.0_lipids_all_201506.itp"

print("Searching for DOPC molecule definition...")
with open(itp_file, 'r') as f:
    for line_num, line in enumerate(f, 1):
        if 'DOPC' in line and line_num < 920:
            print(f"Line {line_num}: '{line.strip()}'")
        if line_num > 920:
            break

print("\nSearching for [bonds] section...")
with open(itp_file, 'r') as f:
    for line_num, line in enumerate(f, 1):
        if '[bonds]' in line and line_num < 950:
            print(f"Line {line_num}: '{line.strip()}'")
        if line_num > 950:
            break

print("\nSearching for bond data...")
with open(itp_file, 'r') as f:
    for line_num, line in enumerate(f, 1):
        if line_num >= 918 and line_num <= 930:
            print(f"Line {line_num}: '{line.strip()}'") 