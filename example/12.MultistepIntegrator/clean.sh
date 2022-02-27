#!/bin/bash

for folder in inputs outputs results
do
    if [ -d $folder ]; then
        rm -rf $folder
    fi
done
