#!/bin/bash

echo "Setting up UPSIDE2 for Apple Silicon Mac (M1/M2/M3/M4)..."
echo `pwd`

# Detect Apple Silicon and set appropriate variables
if [[ $(uname -m) == "arm64" ]]; then
    echo "Detected Apple Silicon (ARM64) architecture"
    ARCH_TYPE="arm64"
else
    echo "Error: This script is designed for Apple Silicon Macs only"
    echo "Detected architecture: $(uname -m)"
    exit 1
fi

# Set the project path dynamically based on current directory
upside_path="$(pwd)"

# Copy and configure the ARM source file
cp source_arm source.sh
sed -i '' "s|UP_PATH|$upside_path|g" source.sh

# Source the environment
source source.sh

# Set up OpenMP environment variables explicitly
LIBOMP_PREFIX=$(brew --prefix libomp)
export OpenMP_CXX_FLAGS="-Xpreprocessor -fopenmp -I$LIBOMP_PREFIX/include"
export OpenMP_CXX_LIB_NAMES="omp"
export OpenMP_C_FLAGS="-Xpreprocessor -fopenmp -I$LIBOMP_PREFIX/include"
export OpenMP_C_LIB_NAMES="omp"
export OpenMP_omp_LIBRARY="$LIBOMP_PREFIX/lib/libomp.dylib"
export CMAKE_PREFIX_PATH="$LIBOMP_PREFIX:$CMAKE_PREFIX_PATH"

echo "OpenMP configuration:"
echo "  LIBOMP_PREFIX: $LIBOMP_PREFIX"
echo "  OpenMP_CXX_FLAGS: $OpenMP_CXX_FLAGS"
echo "  OpenMP_CXX_LIB_NAMES: $OpenMP_CXX_LIB_NAMES"
echo "  OpenMP_omp_LIBRARY: $OpenMP_omp_LIBRARY"

# Clean previous build
echo "Cleaning previous build..."
rm -rf obj/*
cd obj

# Configure with Apple Silicon optimized settings
echo "Configuring CMake for Apple Silicon..."
cmake ../src/ \
  -DEIGEN3_INCLUDE_DIR=$EIGEN_HOME \
  -DCMAKE_OSX_ARCHITECTURES=$ARCH_TYPE \
  -DCMAKE_CXX_FLAGS="-std=c++11" \
  -DCMAKE_C_FLAGS="" \
  -DOpenMP_CXX_FLAGS="$OpenMP_CXX_FLAGS" \
  -DOpenMP_C_FLAGS="$OpenMP_C_FLAGS" \
  -DOpenMP_CXX_LIB_NAMES="$OpenMP_CXX_LIB_NAMES" \
  -DOpenMP_C_LIB_NAMES="$OpenMP_C_LIB_NAMES" \
  -DOpenMP_omp_LIBRARY="$OpenMP_omp_LIBRARY" \
  -DCMAKE_PREFIX_PATH="$CMAKE_PREFIX_PATH"

# Build with verbose output
echo "Building UPSIDE2..."
make VERBOSE=1

# Check if build was successful
if [ $? -eq 0 ]; then
    echo "Build completed successfully!"
    echo "Executables created:"
    ls -la upside*
    echo "Shared library created:"
    ls -la libupside*
else
    echo "Build failed. Please check the error messages above."
    exit 1
fi
 