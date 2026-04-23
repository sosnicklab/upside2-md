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

WORKFLOW_DIR="$(cd "$(dirname "${BASE_WORKFLOW_SCRIPT}")" && pwd)"

autodetect_previous_stage70_file() {
    local outputs_root="${WORKFLOW_DIR}/outputs"
    if [ ! -d "${outputs_root}" ]; then
        return 0
    fi

    python3 - "${outputs_root}" << 'PY'
from pathlib import Path
import sys

root = Path(sys.argv[1])
candidates = []
allowed_names = {"1afo.stage_7.0.up", "1afo.stage_7.0.continue.up"}
for path in root.glob("martini_test_1afo_aabb*/checkpoints/1afo.stage_7.0*.up"):
    if not path.is_file():
        continue
    if path.name not in allowed_names:
        continue
    candidates.append(path)

if not candidates:
    raise SystemExit(0)

candidates.sort(key=lambda path: (path.stat().st_mtime_ns, str(path)), reverse=True)
print(candidates[0])
PY
}

if [ -z "${CONTINUE_STAGE_70_FROM:-}" ]; then
    if [ -n "${PREVIOUS_STAGE7_FILE:-}" ]; then
        export CONTINUE_STAGE_70_FROM="${PREVIOUS_STAGE7_FILE}"
    elif [ -n "${PREVIOUS_RUN_DIR:-}" ]; then
        export CONTINUE_STAGE_70_FROM="${PREVIOUS_RUN_DIR}/checkpoints/1afo.stage_7.0.up"
    elif [ "${AUTO_CONTINUE_FROM_PREVIOUS_RUN:-1}" = "1" ] && [ "${DISABLE_1AFO_AABB_AUTO_CONTINUE:-0}" != "1" ]; then
        auto_continue_file="$(autodetect_previous_stage70_file || true)"
        if [ -n "${auto_continue_file}" ]; then
            export CONTINUE_STAGE_70_FROM="${auto_continue_file}"
        fi
    fi
fi

if [ -n "${CONTINUE_STAGE_70_FROM:-}" ]; then
    export RUN_DIR="${RUN_DIR:-outputs/martini_test_1afo_aabb_continue}"
else
    export RUN_DIR="${RUN_DIR:-outputs/martini_test_1afo_aabb}"
fi

export RUNTIME_PDB_ID="${RUNTIME_PDB_ID:-1afo_aabb}"
export PROTEIN_AA_PDB="${PROTEIN_AA_PDB:-pdb/1AFO.pdb}"

exec "${BASE_WORKFLOW_SCRIPT}" "PDB_ID=${PDB_ID:-1afo}" "$@"
