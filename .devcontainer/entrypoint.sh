#!/bin/bash
set -e

# Source conda setup and activate environment
source /opt/conda/etc/profile.d/conda.sh
conda activate upside2-env

# Set project-specific environment variables for non-interactive scripts.
export UPSIDE_HOME="/upside2-md"
export PATH="$UPSIDE_HOME/py:$UPSIDE_HOME/obj:$PATH"
export PYTHONPATH="$UPSIDE_HOME/py:$PYTHONPATH"

# If no arguments, start an interactive shell
if [ $# -eq 0 ]; then
    exec bash
else
    # Run user's command
    exec "$@"
fi
