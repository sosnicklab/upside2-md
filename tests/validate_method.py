#!/usr/bin/env python

import argparse
import jsonpath
import json
import sys
import subprocess
import os
from typing import List, Tuple

# ANSI color codes for pretty output
class Colors:
    RED = '\033[91m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    MAGENTA = '\033[95m'
    CYAN = '\033[96m'
    WHITE = '\033[97m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'

def print_colored(text: str, color: str = "", bold: bool = False, end: str = "\n"):
    """Print text with color formatting"""
    formatting = ""
    if bold:
        formatting += Colors.BOLD
    if color:
        formatting += color
    
    print(f"{formatting}{text}{Colors.END}", end=end)

def print_header(text: str):
    """Print a formatted header"""
    print()
    print_colored("=" * 60, Colors.CYAN, bold=True)
    print_colored(f" {text} ", Colors.CYAN, bold=True)
    print_colored("=" * 60, Colors.CYAN, bold=True)

def print_test_result(entry_num: int, method: str, status: str, message: str = ""):
    """Print formatted test result"""
    status_colors = {
        "PASS": Colors.GREEN,
        "FAIL": Colors.RED,
        "NEW": Colors.YELLOW,
        "SKIP": Colors.MAGENTA,
        "ERROR": Colors.RED
    }
    
    color = status_colors.get(status, Colors.WHITE)
    print_colored(f"[{status:^5}]", color, bold=True, end=" ")
    print_colored(f"Entry {entry_num:3d} | {method}", Colors.WHITE, end="")
    
    if message:
        print_colored(f" | {message}", Colors.WHITE)
    else:
        print()

def lookup_inputs(data_entry, json_paths):
    results = []
    for json_path in json_paths:
        match = jsonpath.match(f".{json_path}", data_entry)
        if match == None:
            print_colored(f"Warning: No match found for JSON path '{json_path}'", Colors.YELLOW)
        else:
            if hasattr(match, 'obj'):
                results.append(match.obj)
            elif isinstance(match, list) and match:
                results.append(match[0])
            else:
                results.append(match)
    return results

def compare_outputs(expected: str, actual: str) -> Tuple[bool, str]:
    """Compare expected and actual outputs, return (is_match, diff_message)"""
    expected_clean = expected.strip()
    actual_clean = actual.strip()
    
    if expected_clean == actual_clean:
        return True, ""
    
    # Create a simple diff message
    diff_msg = f"\n    Expected: {repr(expected_clean)}\n    Got:      {repr(actual_clean)}"
    return False, diff_msg

def main():
    # Initialize ArgumentParser, disabling the default -h/--help
    parser = argparse.ArgumentParser(description="Run an executable with arguments from a JSON file.", add_help=False)
    
    # Manually add arguments (including our custom --help)
    parser.add_argument('--help', action='help', default=argparse.SUPPRESS,
                        help='Show this help message and exit')
    parser.add_argument('--method', type=str, required=True, help='Module where the function is located')
    parser.add_argument('--args', nargs='*', default=[], help='Positional arguments for the function (JSONPaths)')
    parser.add_argument('--json', type=str, required=True, help='Path to JSON file with keyword arguments')
    parser.add_argument('--force-write', action='store_true', help='Force write output to JSON, overwriting if entry exists')
    
    args = parser.parse_args()

    print_header(f"Testing Method: {args.method}")
    print_colored(f"JSON File: {args.json}", Colors.CYAN)
    print_colored(f"Arguments: {args.args if args.args else 'None'}", Colors.CYAN)
    print_colored(f"Force Write: {'Yes' if args.force_write else 'No'}", Colors.CYAN)
    print()

    try:
        with open(args.json, 'r') as f:
            data = json.load(f)
    except Exception as e:
        print_colored(f"Failed to load JSON file: {e}", Colors.RED, bold=True)
        sys.exit(1)

    data_modified = False
    test_results = {
        'passed': 0,
        'failed': 0,
        'new': 0,
        'errors': 0,
        'skipped': 0,
        'failed_entries': []
    }
    
    total_entries = len(data)
    print_colored(f"Running {total_entries} test entries...", Colors.BLUE, bold=True)
    print()

    for i, entry in enumerate(data):
        test_output_key = f"test_output_{args.method}"
        
        try:
            exe_args_from_json = lookup_inputs(entry, args.args)
            str_exe_args = [str(arg) for arg in exe_args_from_json]

            # Run the subprocess
            executable_path = f"{os.environ.get('UPSIDE_HOME')}/tests/build/testbridge"
            cmd = [executable_path, args.method, *str_exe_args]
            
            output = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            
            if output.returncode != 0:
                print_test_result(i, args.method, "ERROR", f"Exit code {output.returncode}")
                if output.stderr.strip():
                    print_colored(f"    Stderr: {output.stderr.strip()}", Colors.RED)
                test_results['errors'] += 1
                continue

        except FileNotFoundError:
            print_test_result(i, args.method, "ERROR", f"Executable not found: {executable_path}")
            test_results['errors'] += 1
            continue
        except subprocess.TimeoutExpired:
            print_test_result(i, args.method, "ERROR", "Test timed out (30s)")
            test_results['errors'] += 1
            continue
        except Exception as e:
            print_test_result(i, args.method, "ERROR", f"Exception: {str(e)}")
            test_results['errors'] += 1
            continue

        # Handle test output comparison/writing
        if args.force_write:
            if test_output_key not in entry or entry[test_output_key] != output.stdout:
                entry[test_output_key] = output.stdout
                print_test_result(i, args.method, "NEW", "Force-wrote output")
                data_modified = True
                test_results['new'] += 1
            else:
                print_test_result(i, args.method, "PASS", "Force-write: same output")
                entry[test_output_key] = output.stdout
                test_results['passed'] += 1
        else:
            if test_output_key in entry:
                # Compare outputs
                is_match, diff_msg = compare_outputs(entry[test_output_key], output.stdout)
                
                if is_match:
                    print_test_result(i, args.method, "PASS")
                    test_results['passed'] += 1
                else:
                    print_test_result(i, args.method, "FAIL", "Output mismatch")
                    print_colored(diff_msg, Colors.RED)
                    test_results['failed'] += 1
                    test_results['failed_entries'].append(i)
            else:
                # New test case
                entry[test_output_key] = output.stdout
                print_test_result(i, args.method, "NEW", "Added new output")
                data_modified = True
                test_results['new'] += 1
    
    # Save updated JSON if needed
    if data_modified:
        try:
            with open(args.json, 'w') as f:
                json.dump(data, f, indent=2)
            print()
            print_colored(f"Updated {args.json} with test outputs.", Colors.GREEN)
        except Exception as e:
            print_colored(f"Failed to save JSON file: {e}", Colors.RED, bold=True)
            sys.exit(1)

    # Print summary
    print_header("Test Summary")
    
    print_colored(f"Total Entries: {total_entries}", Colors.WHITE, bold=True)
    print_colored(f"Passed:        {test_results['passed']}", Colors.GREEN)
    print_colored(f"Failed:        {test_results['failed']}", Colors.RED)
    print_colored(f"New/Updated:   {test_results['new']}", Colors.YELLOW)
    print_colored(f"Errors:        {test_results['errors']}", Colors.MAGENTA)
    
    if test_results['failed'] > 0:
        print()
        print_colored("Failed Entries:", Colors.RED, bold=True)
        for entry_num in test_results['failed_entries']:
            print_colored(f"  - Entry {entry_num}", Colors.RED)
    
    print()
    
    # Determine overall result
    if test_results['errors'] > 0:
        print_colored("❌ TESTING COMPLETED WITH ERRORS", Colors.RED, bold=True)
        sys.exit(2)  # Different exit code for errors vs failures
    elif test_results['failed'] > 0:
        print_colored("❌ SOME TESTS FAILED", Colors.RED, bold=True)
        sys.exit(1)
    elif test_results['new'] > 0:
        print_colored("✅ ALL TESTS PASSED (with new outputs added)", Colors.YELLOW, bold=True)
    else:
        print_colored("✅ ALL TESTS PASSED", Colors.GREEN, bold=True)

if __name__ == "__main__":
    main()