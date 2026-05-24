#!/bin/bash
#SBATCH --job-name=1afo_hybrid
#SBATCH --output=slurm-%x-%j.out
#SBATCH --time=36:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${UPSIDE_PROJECT_ROOT:-$(cd "${SCRIPT_DIR}/../.." && pwd)}"

if [ -f /etc/profile.d/modules.sh ]; then
    source /etc/profile.d/modules.sh
fi

if command -v module >/dev/null 2>&1; then
    module load python/3.11.9 || true
    module load cmake || true
    module load openmpi || true
    module load "${UPSIDE_HDF5_MODULE:-hdf5/1.14.3}" || true
fi

if [ -f "${PROJECT_ROOT}/.venv/bin/activate" ]; then
    source "${PROJECT_ROOT}/.venv/bin/activate"
fi

export UPSIDE_SKIP_SOURCE_SH="${UPSIDE_SKIP_SOURCE_SH:-1}"
export UPSIDE_HOME="${PROJECT_ROOT}"
export PYTHONUNBUFFERED="${PYTHONUNBUFFERED:-1}"
export PYTHONPATH="${PROJECT_ROOT}/py${PYTHONPATH:+:$PYTHONPATH}"
export PATH="${PROJECT_ROOT}/obj:$PATH"

export PDB_ID="${PDB_ID:-1afo}"
export LIPID_RESOLUTION="${LIPID_RESOLUTION:-coarse}"
export RUNTIME_PDB_ID="${RUNTIME_PDB_ID:-1afo_hybrid}"
export RUN_DIR="${RUN_DIR:-outputs/martini_1afo_hybrid}"
export PROTEIN_AA_PDB="${PROTEIN_AA_PDB:-pdb/1AFO.pdb}"

exec "${SCRIPT_DIR}/run_sim_hybrid.sh" "$@"
