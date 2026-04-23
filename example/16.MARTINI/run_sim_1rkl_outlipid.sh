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
WORKFLOW_DIR="$(cd "$(dirname "${BASE_WORKFLOW_SCRIPT}")" && pwd)"

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

select_latest_stage70_file() {
    local search_root="$1"
    local glob_pattern="$2"
    python3 - "$search_root" "$glob_pattern" << 'PY'
from pathlib import Path
import re
import sys

root = Path(sys.argv[1])
glob_pattern = sys.argv[2]
if not root.exists():
    raise SystemExit(0)

name_pattern = re.compile(r"^1rkl\.stage_7\.(\d+)\.up$")
candidates = []
for path in root.glob(glob_pattern):
    if not path.is_file():
        continue
    match = name_pattern.fullmatch(path.name)
    if not match:
        continue
    candidates.append((int(match.group(1)), path.stat().st_mtime_ns, str(path), path))

if not candidates:
    raise SystemExit(0)

candidates.sort(key=lambda item: (item[0], item[1], item[2]), reverse=True)
print(candidates[0][3])
PY
}

autodetect_previous_stage70_file() {
    select_latest_stage70_file "${WORKFLOW_DIR}/outputs" "martini_test_1rkl_outlipid*/checkpoints/1rkl.stage_7.*.up"
}

resolve_previous_stage70_from_run_dir() {
    local previous_run_dir="$1"
    python3 - "$previous_run_dir" << 'PY'
from pathlib import Path
import re
import sys

root = Path(sys.argv[1]).resolve()
checkpoints = root / "checkpoints"
search_dir = checkpoints if checkpoints.is_dir() else root
name_pattern = re.compile(r"^1rkl\.stage_7\.(\d+)\.up$")
candidates = []
if search_dir.is_dir():
    for path in search_dir.glob("1rkl.stage_7.*.up"):
        if not path.is_file():
            continue
        match = name_pattern.fullmatch(path.name)
        if not match:
            continue
        candidates.append((int(match.group(1)), path.stat().st_mtime_ns, str(path), path))

if not candidates:
    raise SystemExit(0)

candidates.sort(key=lambda item: (item[0], item[1], item[2]), reverse=True)
print(candidates[0][3])
PY
}

derive_run_dir_from_stage70_file() {
    local stage70_file="$1"
    python3 - "$stage70_file" << 'PY'
from pathlib import Path
import sys

path = Path(sys.argv[1]).resolve()
parent = path.parent
if parent.name == "checkpoints":
    print(parent.parent)
else:
    print(parent)
PY
}

# Outside-of-bilayer 1RKL start:
# - place the protein above the upper leaflet,
# - rotate it into a laid-flat orientation,
# - enlarge the box so unfolding has more room than the embedded workflow.
# Continuation options:
# - by default the wrapper auto-detects the newest previous outlipid
#   stage-7 artifact under `outputs/`,
# - set CONTINUE_STAGE_70_FROM directly to a previous stage_7.*.up file, or
# - set PREVIOUS_STAGE7_FILE to that file, or
# - set PREVIOUS_RUN_DIR and the wrapper will use the newest
#   ${PREVIOUS_RUN_DIR}/checkpoints/1rkl.stage_7.*.up
# - set AUTO_CONTINUE_FROM_PREVIOUS_RUN=0 to force a scratch start even
#   when a previous outlipid stage-7 artifact exists.
# Seed options:
# - leave PREP_SEED and SEED unset to let the base workflow generate them
#   randomly per run,
# - set PREP_SEED and/or SEED explicitly for reproducible reruns.
if [ -z "${CONTINUE_STAGE_70_FROM:-}" ]; then
    if [ -n "${PREVIOUS_STAGE7_FILE:-}" ]; then
        export CONTINUE_STAGE_70_FROM="${PREVIOUS_STAGE7_FILE}"
    elif [ -n "${PREVIOUS_RUN_DIR:-}" ]; then
        previous_stage70_file="$(resolve_previous_stage70_from_run_dir "${PREVIOUS_RUN_DIR}" || true)"
        if [ -n "${previous_stage70_file}" ]; then
            export CONTINUE_STAGE_70_FROM="${previous_stage70_file}"
        fi
    elif [ "${AUTO_CONTINUE_FROM_PREVIOUS_RUN:-1}" = "1" ]; then
        auto_continue_file="$(autodetect_previous_stage70_file || true)"
        if [ -n "${auto_continue_file}" ]; then
            export CONTINUE_STAGE_70_FROM="${auto_continue_file}"
        fi
    fi
fi

if [ -n "${CONTINUE_STAGE_70_FROM:-}" ]; then
    if [ -z "${RUN_DIR:-}" ]; then
        derived_run_dir="$(derive_run_dir_from_stage70_file "${CONTINUE_STAGE_70_FROM}" || true)"
        if [ -n "${derived_run_dir}" ]; then
            export RUN_DIR="${derived_run_dir}"
        fi
    fi
fi
export RUN_DIR="${RUN_DIR:-outputs/martini_test_1rkl_outlipid}"

export RUNTIME_PDB_ID="${RUNTIME_PDB_ID:-1rkl_outlipid}"
export BILAYER_PDB="${BILAYER_PDB:-pdb/DOPC.pdb}"
export PROTEIN_PLACEMENT_MODE="${PROTEIN_PLACEMENT_MODE:-outside-top}"
export PROTEIN_ORIENTATION_MODE="${PROTEIN_ORIENTATION_MODE:-lay-flat}"
export PROTEIN_SURFACE_GAP="${PROTEIN_SURFACE_GAP:-6.0}"
export XY_SCALE="${XY_SCALE:-1.35}"
export BOX_PADDING_XY="${BOX_PADDING_XY:-20.0}"
export BOX_PADDING_Z="${BOX_PADDING_Z:-50.0}"
export PROTEIN_LIPID_CUTOFF="${PROTEIN_LIPID_CUTOFF:-5.0}"
export PROTEIN_LIPID_MIN_GAP="${PROTEIN_LIPID_MIN_GAP:-5.0}"
export PROTEIN_LIPID_CUTOFF_MAX="${PROTEIN_LIPID_CUTOFF_MAX:-10.0}"

exec "${BASE_WORKFLOW_SCRIPT}" "$@"
