# Environment for ConDiv training runs.  Source this, do not execute it.
#
# PROJECT_ROOT is the nearest directory above this file that holds py/upside_config.py, so the same
# file works as training/env.sh and as a run directory's copy, unedited. The code and the binary
# always come from this tree: py/ is first on PYTHONPATH, and upside_engine loads the libupside.so
# beside it. The Python is the tree's .venv. On a cluster, load the modules libupside.so and the
# venv's interpreter need (compiler runtime, HDF5, Python) before sourcing this, or add them to the
# run directory's copy.

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && until [ -f py/upside_config.py ] || [ "$PWD" = / ]; do cd ..; done; pwd)"

if [ -f "$PROJECT_ROOT/.venv/bin/activate" ]; then
    source "$PROJECT_ROOT/.venv/bin/activate"
fi

export HDF5_USE_FILE_LOCKING=FALSE
export UPSIDE_HOME="$PROJECT_ROOT"
export PATH="$PROJECT_ROOT/obj:$PATH"
export PYTHONPATH="$PROJECT_ROOT/py${PYTHONPATH:+:$PYTHONPATH}"
export UPSIDE_SKIP_SOURCE_SH=1
