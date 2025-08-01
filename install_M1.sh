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
  -DCMAKE_C_FLAGS=""

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
