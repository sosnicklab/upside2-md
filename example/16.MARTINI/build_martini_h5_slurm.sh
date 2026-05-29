#!/bin/bash
#SBATCH --job-name=martini_h5
#SBATCH --output=slurm-%x-%j.out
#SBATCH --time=36:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G

set -euo pipefail

# Resolve PROJECT_ROOT from the submission directory (Slurm copies the script
# to a spool dir, so BASH_SOURCE[0] is unreliable).  Walk up from
# SLURM_SUBMIT_DIR until we find py/martini_gen_params.py.
if [ -n "${SLURM_SUBMIT_DIR+x}" ]; then
    _cur="${SLURM_SUBMIT_DIR}"
else
    _cur="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
while [ "$_cur" != / ] && [ ! -f "$_cur/py/martini_gen_params.py" ]; do
    _cur="$(dirname "$_cur")"
done
if [ ! -f "$_cur/py/martini_gen_params.py" ]; then
    echo "ERROR: cannot locate project root (no py/martini_gen_params.py found)" >&2
    exit 1
fi
PROJECT_ROOT="$_cur"
unset _cur

if [ -f /etc/profile.d/modules.sh ]; then
    source /etc/profile.d/modules.sh
fi

if command -v module >/dev/null 2>&1; then
    HDF5_MODULE=hdf5/1.14.3
    if [ -n "${HYBRID_SWEEP_HDF5_MODULE+x}" ]; then
        HDF5_MODULE="${HYBRID_SWEEP_HDF5_MODULE}"
    fi
    if [ -n "${UPSIDE_HDF5_MODULE+x}" ]; then
        HDF5_MODULE="${UPSIDE_HDF5_MODULE}"
    fi
    module load python/3.11.9 || true
    module load cmake || true
    module load openmpi || true
    module load "${HDF5_MODULE}" || true
fi

if [ -f "${PROJECT_ROOT}/.venv/bin/activate" ]; then
    source "${PROJECT_ROOT}/.venv/bin/activate"
fi

export UPSIDE_HOME="${PROJECT_ROOT}"
export UPSIDE_SKIP_SOURCE_SH=1
export PYTHONUNBUFFERED=1
export PATH="${PROJECT_ROOT}/obj:$PATH"
export PYTHONPATH="${PROJECT_ROOT}/py${PYTHONPATH:+:$PYTHONPATH}"

if [ -z "${UPSIDE_MARTINI_TABLE_WORKERS+x}" ]; then
    if [ -n "${SLURM_CPUS_PER_TASK+x}" ]; then
        UPSIDE_MARTINI_TABLE_WORKERS="${SLURM_CPUS_PER_TASK}"
    else
        UPSIDE_MARTINI_TABLE_WORKERS=1
    fi
    export UPSIDE_MARTINI_TABLE_WORKERS
fi

# Hidden-bead relaxation is always enabled.  The rigid-geometry table path
# has been removed; every interaction involving DOPC or SC beads is relaxed
# during table construction to produce a physically realistic effective pair
# potential.
UPSIDE_MARTINI_FIT_RELAX_STEPS=50
export UPSIDE_MARTINI_FIT_RELAX_STEPS

echo "Regenerating dry-MARTINI .h5 files under ${PROJECT_ROOT}/parameters/dryMARTINI"
echo "Using ${UPSIDE_MARTINI_TABLE_WORKERS} MARTINI table worker(s), ${UPSIDE_MARTINI_FIT_RELAX_STEPS} fit relax step(s)"

# --fit-relax-steps comes after "$@" so that any caller-supplied value
# cannot accidentally disable relaxation.
python3 "${PROJECT_ROOT}/py/martini_gen_params.py" \
    --upside-home "${PROJECT_ROOT}" \
    --force \
    --workers "${UPSIDE_MARTINI_TABLE_WORKERS}" \
    "$@" \
    --fit-relax-steps "${UPSIDE_MARTINI_FIT_RELAX_STEPS}"
