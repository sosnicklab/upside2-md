#!/bin/bash

upside_path=$(pwd |sed -e 's/\//\\\//g')
echo $upside_path
cp source_sh source.sh
sed -i "s/UP_PATH/$upside_path/g" source.sh

source source.sh

rm -rf obj/*
cd obj

cmake ../src/  -DEIGEN3_INCLUDE_DIR=$EIGEN_HOME
make
