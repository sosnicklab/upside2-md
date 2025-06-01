#!/usr/bin/env python

import subprocess
import re
import os
import json
import traceback

upside_home = os.environ.get("UPSIDE_HOME")

# Instead of function name, use file and line number
source_file = "src/interaction_graph.h"  # Adjust path as needed
line_number = 182  # Replace with actual line number
max_hits = 20

upside_args = [
    "--duration", "1000",
    "--frame-interval", "50", 
    "--temperature", "0.8",
    "--seed", "1",
    f"{upside_home}/example/01.GettingStarted/outputs/simple_test/chig.run.up"
]

def clean_json_values(data):
    """
    Recursively clean JSON data by removing values that are just '{' or other incomplete values.
    """
    if isinstance(data, dict):
        cleaned = {}
        for key, value in data.items():
            if isinstance(value, str):
                # Remove values that are just opening braces or other incomplete markers
                cleaned_value = value.strip()
                if cleaned_value not in ["{", "}", "{}", "(", ")", "()", "[", "]", "[]", "<incomplete type>"]:
                    cleaned[key] = cleaned_value
            elif isinstance(value, (dict, list)):
                cleaned_result = clean_json_values(value)
                if cleaned_result:  # Only add if not empty after cleaning
                    cleaned[key] = cleaned_result
            else:
                cleaned[key] = value
        return cleaned
    elif isinstance(data, list):
        cleaned = []
        for item in data:
            if isinstance(item, str):
                cleaned_item = item.strip()
                if cleaned_item not in ["{", "}", "{}", "(", ")", "()", "[", "]", "[]", "<incomplete type>"]:
                    cleaned.append(cleaned_item)
            elif isinstance(item, (dict, list)):
                cleaned_result = clean_json_values(item)
                if cleaned_result:  # Only add if not empty after cleaning
                    cleaned.append(cleaned_result)
            else:
                cleaned.append(item)
        return cleaned
    else:
        return data

def generate_gdb_script_content(source_file, line_number, log_file_path, max_hits=1):
    """
    Generates GDB script that sets breakpoint by file:line and captures dereferenced values.
    """
    abs_log_file_path = os.path.abspath(log_file_path)

    script_content = f"""set pagination off
set confirm off
set width 0
set print pretty on
set print elements 100
set print max-depth 5
set print null-stop on

set $max_hits = {max_hits}
set $hit_count = 0

set logging file {abs_log_file_path}
set logging redirect on
set logging overwrite on
set logging on

echo Starting breakpoint at {source_file}:{line_number}\\n

break {source_file}:{line_number}
commands
set $hit_count = $hit_count + 1
echo --- GDB_CAPTURE_HIT_START ---\\n
echo HIT_NUMBER: 
print $hit_count
echo LOCATION: Hit at {source_file}:{line_number}\\n

echo === ARGUMENTS ===\\n
info args

echo === DEREFERENCED_ARGUMENTS ===\\n
python
try:
    frame = gdb.selected_frame()
    block = frame.block()
    for symbol in block:
        if symbol.is_argument:
            name = symbol.name
            try:
                value = frame.read_var(name)
                print(f"VARIABLE: {{name}} = {{value}}")
                
                # Try to dereference if it's a pointer
                if str(value.type).endswith('*') or str(value.type).find('*') != -1:
                    try:
                        deref_value = value.dereference()
                        print(f"DEREF: *{{name}} = {{deref_value}}")
                        
                        # For 'this' pointer, try to print member variables
                        if name == 'this':
                            try:
                                obj_type = deref_value.type
                                print(f"THIS_TYPE: {{name}} type: {{obj_type}}")
                                for field in obj_type.fields():
                                    try:
                                        field_value = deref_value[field.name]
                                        print(f"MEMBER: {{name}}->.{{field.name}} = {{field_value}}")
                                        
                                        # Try to dereference member if it's also a pointer
                                        if str(field_value.type).endswith('*') or str(field_value.type).find('*') != -1:
                                            try:
                                                member_deref = field_value.dereference()
                                                print(f"MEMBER_DEREF: *{{name}}->.{{field.name}} = {{member_deref}}")
                                            except:
                                                pass
                                    except Exception as e:
                                        print(f"MEMBER_ERROR: {{name}}->.{{field.name}} = <could not read: {{e}}>")
                            except Exception as e:
                                print(f"THIS_ERROR: Could not read *{{name}} members: {{e}}")
                        
                        # For other pointers, try to access common members if it's a struct/class
                        else:
                            try:
                                obj_type = deref_value.type
                                if hasattr(obj_type, 'fields'):
                                    for field in obj_type.fields()[:5]:  # Limit to first 5 fields
                                        try:
                                            field_value = deref_value[field.name]
                                            print(f"PTR_MEMBER: {{name}}->.{{field.name}} = {{field_value}}")
                                        except:
                                            pass
                            except:
                                pass
                                
                    except Exception as e:
                        print(f"DEREF_ERROR: Could not dereference {{name}}: {{e}}")
                        
            except Exception as e:
                print(f"VAR_ERROR: Could not read {{name}}: {{e}}")
except Exception as e:
    print(f"SCRIPT_ERROR: Error in Python script: {{e}}")
end

echo === LOCAL_VARIABLES ===\\n
info locals

echo === DEREFERENCED_LOCALS ===\\n
python
try:
    frame = gdb.selected_frame()
    block = frame.block()
    for symbol in block:
        if symbol.is_variable and not symbol.is_argument:
            name = symbol.name
            try:
                value = frame.read_var(name)
                print(f"LOCAL: {{name}} = {{value}}")
                
                # Try to dereference if it's a pointer
                if str(value.type).endswith('*') or str(value.type).find('*') != -1:
                    try:
                        deref_value = value.dereference()
                        print(f"LOCAL_DEREF: *{{name}} = {{deref_value}}")
                        
                        # Try to access members if it's a struct/class
                        try:
                            obj_type = deref_value.type
                            if hasattr(obj_type, 'fields'):
                                for field in obj_type.fields()[:5]:  # Limit to first 5 fields
                                    try:
                                        field_value = deref_value[field.name]
                                        print(f"LOCAL_MEMBER: {{name}}->.{{field.name}} = {{field_value}}")
                                    except:
                                        pass
                        except:
                            pass
                            
                    except Exception as e:
                        print(f"LOCAL_DEREF_ERROR: Could not dereference {{name}}: {{e}}")
                        
            except Exception as e:
                print(f"LOCAL_ERROR: Could not read {{name}}: {{e}}")
except Exception as e:
    print(f"LOCAL_SCRIPT_ERROR: Error in Python script: {{e}}")
end

echo === BACKTRACE ===\\n
backtrace 5

echo --- GDB_CAPTURE_HIT_END ---\\n

if $hit_count >= $max_hits
    echo Max hits reached. Stopping execution.\\n
    disable breakpoints
    quit
else
    continue
end
end

run
quit
"""
    return script_content

def find_function_location(executable_path, function_pattern):
    """
    Helper function to find the file:line location of a function using objdump or nm.
    """
    try:
        # Try to find the function using objdump
        result = subprocess.run(['objdump', '-t', executable_path], 
                              capture_output=True, text=True)
        
        # Look for function in symbol table
        for line in result.stdout.split('\n'):
            if function_pattern in line:
                print(f"Found symbol: {line}")
                
        # Alternative: use nm to find symbols
        result = subprocess.run(['nm', '-C', '-l', executable_path], 
                              capture_output=True, text=True)
        
        for line in result.stdout.split('\n'):
            if function_pattern in line and ':' in line:
                parts = line.split()
                if len(parts) >= 3:
                    location = parts[-1]  # Usually the last part is file:line
                    if ':' in location:
                        file_part, line_part = location.rsplit(':', 1)
                        try:
                            line_num = int(line_part)
                            print(f"Found function at: {file_part}:{line_num}")
                            return file_part, line_num
                        except ValueError:
                            continue
                            
    except Exception as e:
        print(f"Error finding function location: {e}")
    
    return None, None

def list_source_files(executable_path):
    """
    List source files that were compiled into the executable.
    """
    try:
        result = subprocess.run(['objdump', '-g', executable_path], 
                              capture_output=True, text=True)
        
        files = set()
        for line in result.stdout.split('\n'):
            if 'DW_AT_name' in line and ('.cpp' in line or '.h' in line):
                # Extract filename from debug info
                match = re.search(r'DW_AT_name\s*:\s*(.+\.(cpp|h|hpp|cc))', line)
                if match:
                    files.add(match.group(1))
        
        return sorted(files)
    except Exception as e:
        print(f"Error listing source files: {e}")
        return []

def parse_gdb_log(log_file_path):
    """
    Parse the GDB log file to extract captured data with improved handling of all output types.
    """
    all_hits_data = []
    current_hit_data = None
    current_section = None
    
    try:
        with open(log_file_path, "r") as f:
            content = f.read()
            
        # Split into hits
        hits = content.split("--- GDB_CAPTURE_HIT_START ---")[1:]
        
        for hit_content in hits:
            if "--- GDB_CAPTURE_HIT_END ---" not in hit_content:
                continue
                
            hit_content = hit_content.split("--- GDB_CAPTURE_HIT_END ---")[0]
            
            hit_data = {
                "hit_number": None,
                "location": None,
                "args": {},
                "locals": {},
                "backtrace": [],
                "dereferenced_args": {},
                "dereferenced_locals": {},
                "this_members": {},
                "pointer_members": {},
                "errors": []
            }
            
            lines = hit_content.split('\n')
            i = 0
            
            while i < len(lines):
                line = lines[i].strip()
                
                # Parse hit number
                if line.startswith("HIT_NUMBER:"):
                    i += 1
                    if i < len(lines) and lines[i].strip().startswith("$"):
                        hit_data["hit_number"] = lines[i].strip().split("=")[-1].strip()
                    elif i < len(lines) and lines[i].strip().isdigit():
                        hit_data["hit_number"] = lines[i].strip()
                    i += 1
                    continue
                    
                # Parse location
                if line.startswith("LOCATION:"):
                    hit_data["location"] = line.replace("LOCATION:", "").strip()
                    i += 1
                    continue
                
                # Parse sections
                if line == "=== ARGUMENTS ===":
                    current_section = "args"
                    i += 1
                    # Parse regular arguments
                    while i < len(lines) and not lines[i].strip().startswith("==="):
                        arg_line = lines[i].strip()
                        if arg_line and arg_line not in ["No arguments."]:
                            # Handle argument lines like "arg_name = arg_value"
                            if "=" in arg_line:
                                name, value = arg_line.split("=", 1)
                                hit_data["args"][name.strip()] = value.strip()
                        i += 1
                    continue
                    
                elif line == "=== DEREFERENCED_ARGUMENTS ===":
                    current_section = "dereferenced_args"
                    i += 1
                    # Parse dereferenced arguments with prefixes
                    while i < len(lines) and not lines[i].strip().startswith("==="):
                        arg_line = lines[i].strip()
                        if arg_line:
                            if arg_line.startswith("VARIABLE:"):
                                # Regular variable
                                content = arg_line.replace("VARIABLE:", "").strip()
                                if "=" in content:
                                    name, value = content.split("=", 1)
                                    hit_data["dereferenced_args"][name.strip()] = value.strip()
                            elif arg_line.startswith("DEREF:"):
                                # Dereferenced pointer
                                content = arg_line.replace("DEREF:", "").strip()
                                if "=" in content:
                                    name, value = content.split("=", 1)
                                    hit_data["dereferenced_args"][name.strip()] = value.strip()
                            elif arg_line.startswith("THIS_TYPE:"):
                                # This pointer type info
                                content = arg_line.replace("THIS_TYPE:", "").strip()
                                hit_data["dereferenced_args"]["this_type"] = content
                            elif arg_line.startswith("MEMBER:"):
                                # This pointer members
                                content = arg_line.replace("MEMBER:", "").strip()
                                if "=" in content:
                                    name, value = content.split("=", 1)
                                    hit_data["this_members"][name.strip()] = value.strip()
                            elif arg_line.startswith("MEMBER_DEREF:"):
                                # Dereferenced this members
                                content = arg_line.replace("MEMBER_DEREF:", "").strip()
                                if "=" in content:
                                    name, value = content.split("=", 1)
                                    hit_data["this_members"][name.strip() + "_deref"] = value.strip()
                            elif arg_line.startswith("PTR_MEMBER:"):
                                # Other pointer members
                                content = arg_line.replace("PTR_MEMBER:", "").strip()
                                if "=" in content:
                                    name, value = content.split("=", 1)
                                    hit_data["pointer_members"][name.strip()] = value.strip()
                            elif arg_line.startswith(("MEMBER_ERROR:", "DEREF_ERROR:", "VAR_ERROR:", "THIS_ERROR:", "SCRIPT_ERROR:")):
                                # Errors
                                hit_data["errors"].append(arg_line)
                        i += 1
                    continue
                    
                elif line == "=== LOCAL_VARIABLES ===":
                    current_section = "locals"
                    i += 1
                    # Parse regular local variables
                    while i < len(lines) and not lines[i].strip().startswith("==="):
                        local_line = lines[i].strip()
                        if local_line and local_line not in ["No locals."]:
                            if "=" in local_line:
                                name, value = local_line.split("=", 1)
                                hit_data["locals"][name.strip()] = value.strip()
                        i += 1
                    continue
                    
                elif line == "=== DEREFERENCED_LOCALS ===":
                    current_section = "dereferenced_locals"
                    i += 1
                    # Parse dereferenced local variables
                    while i < len(lines) and not lines[i].strip().startswith("==="):
                        local_line = lines[i].strip()
                        if local_line:
                            if local_line.startswith("LOCAL:"):
                                content = local_line.replace("LOCAL:", "").strip()
                                if "=" in content:
                                    name, value = content.split("=", 1)
                                    hit_data["dereferenced_locals"][name.strip()] = value.strip()
                            elif local_line.startswith("LOCAL_DEREF:"):
                                content = local_line.replace("LOCAL_DEREF:", "").strip()
                                if "=" in content:
                                    name, value = content.split("=", 1)
                                    hit_data["dereferenced_locals"][name.strip()] = value.strip()
                            elif local_line.startswith("LOCAL_MEMBER:"):
                                content = local_line.replace("LOCAL_MEMBER:", "").strip()
                                if "=" in content:
                                    name, value = content.split("=", 1)
                                    hit_data["pointer_members"][name.strip()] = value.strip()
                            elif local_line.startswith(("LOCAL_DEREF_ERROR:", "LOCAL_ERROR:", "LOCAL_SCRIPT_ERROR:")):
                                hit_data["errors"].append(local_line)
                        i += 1
                    continue
                    
                elif line == "=== BACKTRACE ===":
                    current_section = "backtrace"
                    i += 1
                    # Parse backtrace
                    while i < len(lines) and not lines[i].strip().startswith("---"):
                        bt_line = lines[i].strip()
                        if bt_line and bt_line.startswith("#"):
                            hit_data["backtrace"].append(bt_line)
                        i += 1
                    continue
                
                i += 1
            
            all_hits_data.append(hit_data)

    except FileNotFoundError:
        print(f"Error: GDB log file not found: {log_file_path}")
        return []
    except Exception as e:
        print(f"Error parsing GDB log file: {e}")
        traceback.print_exc()
        return []
        
    return all_hits_data

def capture_function_inputs_with_gdb(
    executable_path,
    source_file,
    line_number,
    program_args=None,
    max_hits=1,
    gdb_temp_script_path=f"{upside_home}/tests/records/temp_capture_args.gdb",
    gdb_log_file="gdb_function_args.log"
):
    """
    Captures data at a specific file:line location using GDB.
    """
    if not os.path.exists(executable_path):
        print(f"Error: Executable not found at '{executable_path}'")
        return None
    
    log_dir = os.path.dirname(os.path.abspath(gdb_log_file))
    os.makedirs(log_dir, exist_ok=True)

    gdb_script_content = generate_gdb_script_content(source_file, line_number, gdb_log_file, max_hits)
    
    with open(gdb_temp_script_path, "w") as f:
        f.write(gdb_script_content)

    if program_args:
        gdb_command = [
            "gdb",
            "--batch",
            "-nx", 
            "-x", gdb_temp_script_path,
            "--args",
            executable_path
        ] + program_args
    else:
        gdb_command = [
            "gdb",
            "--batch",
            "-nx",
            "-x", gdb_temp_script_path,
            executable_path
        ]

    print(f"Running GDB with command: {' '.join(gdb_command)}")
    
    try:
        timeout_seconds = 120  # Increased timeout
        process = subprocess.run(
            gdb_command,
            capture_output=True, 
            text=True, 
            check=False,
            timeout=timeout_seconds
        )
        
        print("GDB Process STDOUT:\n", process.stdout)
        if process.stderr:
            print("GDB Process STDERR:\n", process.stderr)
        if process.returncode != 0:
            print(f"Warning: GDB process exited with code {process.returncode}")

    except subprocess.TimeoutExpired:
        print(f"GDB command timed out after {timeout_seconds} seconds.")
        return None
    except Exception as e:
        print(f"An error occurred while running GDB: {e}")
        return None

    captured_data = parse_gdb_log(gdb_log_file)
    return captured_data

def save_to_json(data, filename):
    """Save captured data to JSON file after cleaning."""
    try:
        # Clean the data before saving
        cleaned_data = clean_json_values(data)
        
        with open(filename, 'w') as f:
            json.dump(cleaned_data, f, indent=2, default=str)
        print(f"Cleaned data saved to {filename}")
    except Exception as e:
        print(f"Error saving to JSON: {e}")

if __name__ == "__main__":
    cpp_executable = f"{upside_home}/obj/upside"
    function_pattern = "PairlistComputation"  # Pattern to search for
    
    if not os.path.exists(cpp_executable):
        print(f"Executable '{cpp_executable}' not found.")
        exit(1)
    
    # # Try to automatically find the function location
    # print("Searching for function location...")
    # found_file, found_line = find_function_location(cpp_executable, function_pattern)
    
    # if found_file and found_line:
    #     print(f"Found function at: {found_file}:{found_line}")
    #     source_file = found_file
    #     line_number = found_line
    # else:
    #     print("Could not automatically find function. Using manual specification.")
    #     print("Available source files in executable:")
    #     files = list_source_files(cpp_executable)
    #     for f in files[:10]:  # Show first 10
    #         print(f"  {f}")
        
        # Manual specification - you'll need to find the right file and line
        source_file = source_file  # UPDATE THIS
        line_number = line_number  # UPDATE THIS
        print(f"Using manual location: {source_file}:{line_number}")

    captured_inputs = capture_function_inputs_with_gdb(
        executable_path=cpp_executable,
        source_file=source_file,
        line_number=line_number,
        program_args=upside_args,
        max_hits=max_hits,
        gdb_log_file=f"{upside_home}/tests/records/line_{line_number}_inputs.log"
    )

    if captured_inputs:
        print("\n--- Successfully Captured Data ---")
        
        # Save to JSON with cleaning
        json_filename = f"line_{line_number}_capture.json"
        save_to_json(captured_inputs, f'{upside_home}/tests/records/{json_filename}')
        
        # Print summary
        for i, hit_data in enumerate(captured_inputs):
            hit_num = hit_data.get('hit_number', 'Unknown')
            print(f"Hit {hit_num}:")
            print(f"  Location: {hit_data.get('location', 'Unknown')}")
            print(f"  Arguments: {len(hit_data['args'])} items")
            print(f"  Dereferenced Args: {len(hit_data['dereferenced_args'])} items")
            print(f"  'this' Members: {len(hit_data['this_members'])} items")
            print(f"  Pointer Members: {len(hit_data['pointer_members'])} items")
            print(f"  Local Variables: {len(hit_data['locals'])} items")
            print(f"  Dereferenced Locals: {len(hit_data['dereferenced_locals'])} items")
            print(f"  Backtrace: {len(hit_data['backtrace'])} frames")
            print(f"  Errors: {len(hit_data['errors'])} items")
            
            if hit_data['errors']:
                print("  Error details:")
                for error in hit_data['errors'][:3]:  # Show first 3 errors
                    print(f"    {error}")
            print("-" * 40)
            
        print(f"\nCleaned data saved to {json_filename}")
    else:
        print("No data captured.")