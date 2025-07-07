# Testing Framework

This directory provides a minimal, general-purpose framework for testing C++ classes and methods with real runtime data.

## Features

- **Minimal source changes:** To test protected/private members, only a `friend` declaration is needed in your class.
- **General:** Not limited to any specific class. Add new test methods in `testbridge.cpp` for any class.
- **GDB-based input capture:** `gdb_input_capture.py` can extract arguments, locals, and dereferenced arrays at any code location. You can tag custom GDB commands to capture specific arrays or structures.
- **Automated validation:** `validate_method.py` runs your test executable with JSON-captured inputs and compares outputs.
- **Easy extension:** Add new test methods, capture new data, or add new test cases with minimal effort.

## How to Use

1. **Build the test executable:**
    ```bash
    ./build.sh
    ```
2. **Run all tests:**
    ```bash
    ./test_runner.sh
    ```

## Adding/Expanding Tests

1. **Expose internals:** Add a `friend` declaration for your test bridge in the class you want to test.
2. **Add test logic:** Implement new test methods in `testbridge.cpp`.
3. **Capture real inputs:** Use `gdb_input_capture.py` to generate JSON test data from a running executable. You can add custom GDB commands to extract arrays or other data.
4. **Create/extend test cases:** Place JSON files in `tests/records/`.
5. **Validate:** Run `validate_method.py` with your method and JSON file. Use `--force-write` to update outputs.

## Minimal Workflow

- Only add a `friend` declaration to access protected members.
- Add/modify test logic in `testbridge.cpp`.
- Use GDB capture for real data.
- Validate with `validate_method.py`.
