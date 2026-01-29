
# Style and Conventions

## C++ Code Style
- Uses C++11 standard
- Header files have .h extension
- Implementation files have .cpp extension
- Classes and structs use PascalCase (e.g., DerivEngine)
- Methods and functions use snake_case (e.g., compute_rotamer_pos)
- Variables use snake_case (e.g., frame_interval)
- Constants use uppercase with underscores (e.g., MAX_FRAMES)
- Use of Eigen library for matrix operations
- OpenMP for parallelization

## Python Code Style
- Uses Python 3.7+
- Scripts have .py extension
- Functions and variables use snake_case (e.g., load_upside_traj)
- Classes use PascalCase (e.g., UpsideEngine)
- Modules and packages use lowercase with underscores
- Docstrings follow NumPy/SciPy style
- Use of standard libraries and common scientific computing packages

## Build Conventions
- CMake for cross-platform build system
- Out-of-source build in obj/ directory
- Installation scripts for different platforms (install.sh, install_M1.sh)
- Virtual environment for Python dependencies (.venv)

## Configuration Files
- Simulation configurations stored in HDF5 (.up) files
- Parameters in HDF5 (.h5) files and text files
- Environment variables set in source.sh
