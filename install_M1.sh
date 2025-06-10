#!/bin/bash

echo `pwd`
upside_path="/Users/yinhan/Documents/upside2-md"
cp source_arm source.sh
#sed -i '' "s/UP_PATH/$upside_path/g" source.sh
sed -i '' "s|UP_PATH|$upside_path|g" source.sh

source source.sh

rm -rf obj/*
cd obj

cmake ../src/  -DEIGEN3_INCLUDE_DIR=$EIGEN_HOME
make VERBOSE=1
