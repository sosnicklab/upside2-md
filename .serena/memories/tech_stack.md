
# Tech Stack

## C++ Core (src/)
- **Language**: C++11
- **Key Libraries**:
  - Eigen 3.0+ (matrix computations)
  - HDF5 1.8+ (file I/O)
  - OpenMP (parallelism)
  - TCLAP (command line parsing)

## Python Layer (py/)
- **Language**: Python 3.7+
- **Key Libraries**:
  - NumPy
  - SciPy
  - PyTables (HDF5 support)
  - ProDy (PDB parsing)
  - MDTraj (trajectory analysis)
  - Pandas
  - H5py
  - Pymbar (MBAR reweighting)

## Build System
- CMake 2.8+
- Apple Clang (for M1/M2/M3/M4 Macs)

## File Formats
- HDF5 (.up files) for simulation configurations and trajectories
- PDB for initial structures
- FASTA for sequences
- VTF for visualization in VMD
