#!/bin/bash
#SBATCH --job-name=1rkl_outlipid
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
        "${UPSIDE_PROJECT_ROOT:-}" \
        "${SLURM_SUBMIT_DIR:-}" \
        "${SLURM_SUBMIT_DIR:-}/example/16.MARTINI" \
        "${PWD}" \
        "${PWD}/example/16.MARTINI" \
        "${SCRIPT_DIR}"
    do
        if [ -n "${candidate_dir}" ] && [ -f "${candidate_dir}/run_sim_1rkl.sh" ]; then
            BASE_WORKFLOW_SCRIPT="${candidate_dir}/run_sim_1rkl.sh"
            break
        fi
    done
fi

if [ -z "${BASE_WORKFLOW_SCRIPT}" ]; then
    echo "ERROR: could not locate run_sim_1rkl.sh." >&2
    echo "Set BASE_WORKFLOW_SCRIPT explicitly or submit from the repo root or example/16.MARTINI." >&2
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

# Outside-of-bilayer 1RKL start:
# - place the protein above the upper leaflet,
# - rotate it into a laid-flat orientation,
# - enlarge the box so unfolding has more room than the embedded workflow.
# Continuation options:
# - set CONTINUE_STAGE_70_FROM directly to a previous stage_7.0.up file, or
# - set PREVIOUS_STAGE7_FILE to that file, or
# - set PREVIOUS_RUN_DIR and the wrapper will use
#   ${PREVIOUS_RUN_DIR}/checkpoints/1rkl.stage_7.0.up
# Seed options:
# - leave PREP_SEED and SEED unset to let the base workflow generate them
#   randomly per run,
# - set PREP_SEED and/or SEED explicitly for reproducible reruns.
if [ -z "${CONTINUE_STAGE_70_FROM:-}" ]; then
    if [ -n "${PREVIOUS_STAGE7_FILE:-}" ]; then
        export CONTINUE_STAGE_70_FROM="${PREVIOUS_STAGE7_FILE}"
    elif [ -n "${PREVIOUS_RUN_DIR:-}" ]; then
        export CONTINUE_STAGE_70_FROM="${PREVIOUS_RUN_DIR}/checkpoints/1rkl.stage_7.0.up"
    fi
fi

if [ -n "${CONTINUE_STAGE_70_FROM:-}" ]; then
    export RUN_DIR="${RUN_DIR:-outputs/martini_test_1rkl_outlipid_continue}"
else
    export RUN_DIR="${RUN_DIR:-outputs/martini_test_1rkl_outlipid}"
fi

export RUNTIME_PDB_ID="${RUNTIME_PDB_ID:-1rkl_outlipid}"
export PROTEIN_PLACEMENT_MODE="${PROTEIN_PLACEMENT_MODE:-outside-top}"
export PROTEIN_ORIENTATION_MODE="${PROTEIN_ORIENTATION_MODE:-lay-flat}"
export PROTEIN_SURFACE_GAP="${PROTEIN_SURFACE_GAP:-6.0}"
export XY_SCALE="${XY_SCALE:-1.35}"
export BOX_PADDING_XY="${BOX_PADDING_XY:-20.0}"
export BOX_PADDING_Z="${BOX_PADDING_Z:-50.0}"
export PROTEIN_LIPID_CUTOFF="${PROTEIN_LIPID_CUTOFF:-5.0}"
export PROTEIN_LIPID_MIN_GAP="${PROTEIN_LIPID_MIN_GAP:-5.0}"
export PROTEIN_LIPID_CUTOFF_MAX="${PROTEIN_LIPID_CUTOFF_MAX:-10.0}"
if [ -n "${CONTINUE_STAGE_70_FROM:-}" ]; then
    export CONTINUE_STAGE_70_OUTPUT="${CONTINUE_STAGE_70_OUTPUT:-${RUN_DIR}/checkpoints/1rkl.stage_7.0.continue.up}"
fi

exec "${BASE_WORKFLOW_SCRIPT}" "$@"
