#!/bin/bash
mkdir -p build
rm -rf build/*
cd build
cmake .. -DEIGEN3_INCLUDE_DIR=$EIGEN_HOME
make

