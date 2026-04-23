#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_WORKFLOW_SCRIPT="${BASE_WORKFLOW_SCRIPT:-}"

if [ -z "${BASE_WORKFLOW_SCRIPT}" ]; then
    for candidate_dir in \
        "${UPSIDE_PROJECT_ROOT:-}/example/16.MARTINI" \
        "${UPSIDE_PROJECT_ROOT:-}" \
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
    echo "Set BASE_WORKFLOW_SCRIPT explicitly or run from the repo root or example/16.MARTINI." >&2
    exit 1
fi

export RUN_DIR="${RUN_DIR:-outputs/martini_test_1afo_aabb}"
export RUNTIME_PDB_ID="${RUNTIME_PDB_ID:-1afo_aabb}"
export PROTEIN_AA_PDB="${PROTEIN_AA_PDB:-pdb/1AFO.pdb}"

exec "${BASE_WORKFLOW_SCRIPT}" "PDB_ID=${PDB_ID:-1afo}" "$@"
