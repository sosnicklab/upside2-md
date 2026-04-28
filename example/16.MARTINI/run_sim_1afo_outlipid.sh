#!/bin/bash
#SBATCH --job-name=1afo_outlipid
#SBATCH --output=slurm-%x-%j.out
#SBATCH --time=36:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_WORKFLOW_SCRIPT="${BASE_WORKFLOW_SCRIPT:-}"

if [ -z "${BASE_WORKFLOW_SCRIPT}" ]; then
    for candidate_dir in \
        "${UPSIDE_PROJECT_ROOT:-}/example/16.MARTINI" \
        "${SLURM_SUBMIT_DIR:-}" \
        "${SLURM_SUBMIT_DIR:-}/example/16.MARTINI" \
        "${PWD}" \
        "${PWD}/example/16.MARTINI" \
        "${SCRIPT_DIR}"
    do
        if [ -n "${candidate_dir}" ] && [ -f "${candidate_dir}/run_sim_1afo.sh" ]; then
            BASE_WORKFLOW_SCRIPT="${candidate_dir}/run_sim_1afo.sh"
            break
        fi
    done
fi

if [ -z "${BASE_WORKFLOW_SCRIPT}" ]; then
    echo "ERROR: could not locate run_sim_1afo.sh." >&2
    exit 1
fi

PROJECT_ROOT="${UPSIDE_PROJECT_ROOT:-$(cd "$(dirname "${BASE_WORKFLOW_SCRIPT}")/../.." && pwd)}"

if [ -f /etc/profile.d/modules.sh ]; then
    source /etc/profile.d/modules.sh
fi

if command -v module >/dev/null 2>&1; then
    module load python/3.11.9 || true
    module load cmake || true
    module load openmpi || true
    module load "${UPSIDE_HDF5_MODULE:-${HYBRID_SWEEP_HDF5_MODULE:-hdf5/1.14.3}}" || true
fi

if [ -f "${PROJECT_ROOT}/.venv/bin/activate" ]; then
    source "${PROJECT_ROOT}/.venv/bin/activate"
fi

export UPSIDE_SKIP_SOURCE_SH="${UPSIDE_SKIP_SOURCE_SH:-1}"
export UPSIDE_HOME="${PROJECT_ROOT}"
export PYTHONUNBUFFERED="${PYTHONUNBUFFERED:-1}"
export PYTHONPATH="${PROJECT_ROOT}/py${PYTHONPATH:+:$PYTHONPATH}"
export PATH="${PROJECT_ROOT}/obj:$PATH"

export RUN_DIR="${RUN_DIR:-outputs/martini_test_1afo_outlipid}"
export RUNTIME_PDB_ID="${RUNTIME_PDB_ID:-1afo_outlipid}"
export PROTEIN_PLACEMENT_MODE="${PROTEIN_PLACEMENT_MODE:-outside-top}"
export PROTEIN_ORIENTATION_MODE="${PROTEIN_ORIENTATION_MODE:-lay-flat}"
export PROTEIN_SURFACE_GAP="${PROTEIN_SURFACE_GAP:-6.0}"
export XY_SCALE="${XY_SCALE:-1.35}"
export BOX_PADDING_XY="${BOX_PADDING_XY:-20.0}"
export BOX_PADDING_Z="${BOX_PADDING_Z:-50.0}"
export PROTEIN_LIPID_CUTOFF="${PROTEIN_LIPID_CUTOFF:-5.0}"
export PROTEIN_LIPID_MIN_GAP="${PROTEIN_LIPID_MIN_GAP:-5.0}"
export PROTEIN_LIPID_CUTOFF_MAX="${PROTEIN_LIPID_CUTOFF_MAX:-10.0}"
export AUTO_CONTINUE_FROM_PREVIOUS_RUN="${AUTO_CONTINUE_FROM_PREVIOUS_RUN:-1}"
export AUTO_CONTINUE_GLOB="${AUTO_CONTINUE_GLOB:-martini_test_1afo_outlipid*/checkpoints/1afo.stage_7*.up}"

exec "${BASE_WORKFLOW_SCRIPT}" "$@"
