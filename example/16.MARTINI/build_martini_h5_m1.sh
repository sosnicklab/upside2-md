#!/bin/bash
set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

if [ -f "${PROJECT_ROOT}/.venv/bin/activate" ]; then
    source "${PROJECT_ROOT}/.venv/bin/activate"
else
    echo "ERROR: missing Python environment: ${PROJECT_ROOT}/.venv" >&2
    exit 1
fi

source "${PROJECT_ROOT}/source.sh"

set -u

export UPSIDE_HOME="${PROJECT_ROOT}"
export PYTHONUNBUFFERED=1
export PYTHONPATH="${PROJECT_ROOT}/py${PYTHONPATH:+:$PYTHONPATH}"
export PATH="${PROJECT_ROOT}/obj:$PATH"

if [ -z "${UPSIDE_MARTINI_TABLE_WORKERS+x}" ]; then
    if command -v sysctl >/dev/null 2>&1; then
        UPSIDE_MARTINI_TABLE_WORKERS="$(sysctl -n hw.logicalcpu)"
    else
        UPSIDE_MARTINI_TABLE_WORKERS="$(getconf _NPROCESSORS_ONLN)"
    fi
    export UPSIDE_MARTINI_TABLE_WORKERS
fi

echo "Regenerating dry-MARTINI .h5 files under ${PROJECT_ROOT}/parameters/dryMARTINI"
echo "Using ${UPSIDE_MARTINI_TABLE_WORKERS} MARTINI table worker(s)"

python3 "${PROJECT_ROOT}/py/martini_gen_params.py" \
    --upside-home "${PROJECT_ROOT}" \
    --force \
    --workers "${UPSIDE_MARTINI_TABLE_WORKERS}" \
    "$@"
