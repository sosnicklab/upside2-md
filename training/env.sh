# Environment for ConDiv training runs.  Source this, do not execute it.
#
# PROJECT_ROOT is the nearest directory above this file that holds py/upside_config.py, so the same
# file works as training/env.sh and as a run directory's copy, unedited. The code and the binary
# always come from this tree: py/ is first on PYTHONPATH, and upside_engine loads the libupside.so
# beside it. Only the Python differs by cluster, with the same package versions on both (numpy
# 1.23.5, scipy 1.13.1, tables 3.8.0, torch 2.6.0+cpu), so a checkpoint written on one cluster
# loads on the other:
#
#   midway2  this tree's .venv, built from midway2's python/3.9.18 module. gcc/10.1.0 is required
#            (the system libstdc++ lacks GLIBCXX_3.4.20 and libupside.so will not load without
#            it), and the python module must be loaded BEFORE the venv is activated, or the
#            interpreter cannot find libpython3.9.so.1.0.
#   midway3  the shared /beagle3 deployment's (env_shared.sh): an el7 interpreter and venv on
#            shared storage. This tree's .venv points at midway2's /software, which midway3 lacks.
#   locally  the repo's .venv.

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && until [ -f py/upside_config.py ] || [ "$PWD" = / ]; do cd ..; done; pwd)"
SHARED_ENV=/beagle3/trsosnic/yinhan/upside2-md/env_shared.sh

# The cluster is told by whether this tree's venv interpreter exists here, which is exactly what
# matters: /software/modules/init/bash exists on midway3 too, so it does not identify midway2.
if [ -f /software/modules/init/bash ] && [ -x "$(readlink -f "$PROJECT_ROOT/.venv/bin/python3")" ]; then
    source /software/modules/init/bash
    module load gcc/10.1.0
    module load python/3.9.18 hdf5/1.14.3+oneapi-2023.1
    source "$PROJECT_ROOT/.venv/bin/activate"
elif [ -f "$SHARED_ENV" ]; then
    source "$SHARED_ENV"
elif [ -f "$PROJECT_ROOT/.venv/bin/activate" ]; then
    source "$PROJECT_ROOT/.venv/bin/activate"
fi

export HDF5_USE_FILE_LOCKING=FALSE
export UPSIDE_HOME="$PROJECT_ROOT"
export PATH="$PROJECT_ROOT/obj:$PATH"
# $PROJECT_ROOT/training carries the trainer's own helpers (rama_basin); py/ is shared Upside
# infrastructure and training-specific code does not belong there.
export PYTHONPATH="$PROJECT_ROOT/py:$PROJECT_ROOT/training${PYTHONPATH:+:$PYTHONPATH}"
export UPSIDE_SKIP_SOURCE_SH=1
