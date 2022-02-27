#!/bin/bash
module load cmake
export UPSIDE_HOME=/home/pengxd/data/source/upside-md
#================== upside ======================
export PATH=$UPSIDE_HOME/obj:$PATH
export PYTHONPATH=$UPSIDE_HOME/py:$PYTHONPATH

#================== python ======================
export MY_PYTHON=/project2/trsosnic/software/Anaconda3
export PATH="$MY_PYTHON/bin:$PATH"

#================== eigen ======================
export EIGEN_HOME=/project2/trsosnic/software/Eigen-3.3.5
