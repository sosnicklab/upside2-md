# Environment for ConDiv training runs.  Source this, do not execute it.
#
# PROJECT_ROOT is derived from this file's location, so a training directory can be copied
# anywhere inside the repo without editing paths.
#
# midway2 notes, both learned the hard way:
#   * gcc/10.1.0 is required -- the system libstdc++ lacks GLIBCXX_3.4.20 and libupside.so
#     will not load without it.
#   * load the python module BEFORE activating .venv, or the venv interpreter cannot find
#     libpython3.9.so.1.0.

TRAIN_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$TRAIN_DIR/.." && pwd)"

if [ -f /software/modules/init/bash ]; then
    source /software/modules/init/bash
    module load gcc/10.1.0
    module load python/3.9.18 hdf5/1.14.3+oneapi-2023.1
elif [ -f /etc/profile.d/modules.sh ]; then
    source /etc/profile.d/modules.sh
    module load python/3.11.9 || true
    module load hdf5/1.14.3   || true
fi

export HDF5_USE_FILE_LOCKING=FALSE

if [ -f "$PROJECT_ROOT/.venv/bin/activate" ]; then
    source "$PROJECT_ROOT/.venv/bin/activate"
fi

export UPSIDE_HOME="$PROJECT_ROOT"
export PATH="$PROJECT_ROOT/obj:$PATH"
export PYTHONPATH="$PROJECT_ROOT/py${PYTHONPATH:+:$PYTHONPATH}"
export UPSIDE_SKIP_SOURCE_SH=1
