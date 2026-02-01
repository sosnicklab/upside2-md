#!/bin/bash

echo "Setting up UPSIDE2 for Apple Silicon Mac (M1/M2/M3/M4) with Apple Clang..."
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

# Remove existing source.sh and generate a fresh one from template
echo "Generating fresh source.sh from template..."
if [ -f "source.sh" ]; then
    echo "Backing up existing source.sh to source.sh.backup"
    cp source.sh source.sh.backup
fi
rm -f source.sh
cp source_arm source.sh
sed -i '' "s|UP_PATH|$upside_path|g" source.sh
echo "Generated source.sh with path: $upside_path"

# Source the environment
source source.sh

# Clean previous build
echo "Cleaning previous build..."
rm -rf obj/*
cd obj

# Configure with Apple Silicon optimized settings using Apple Clang
echo "Configuring CMake for Apple Silicon with Apple Clang..."
cmake ../src/ \
  -DEIGEN3_INCLUDE_DIR=$EIGEN_HOME \
  -DCMAKE_OSX_ARCHITECTURES=$ARCH_TYPE \
  -DCMAKE_CXX_COMPILER=/usr/bin/g++ \
  -DCMAKE_C_COMPILER=/usr/bin/gcc \
  -DCMAKE_CXX_FLAGS="-std=c++11" \
  -DCMAKE_C_FLAGS="" \
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
