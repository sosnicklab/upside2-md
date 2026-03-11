#!/bin/bash
set -euo pipefail

resolve_paths() {
  local candidate

  if [ -n "${CONDIV_PROJECT_ROOT:-}" ] && [ -f "${CONDIV_PROJECT_ROOT}/source.sh" ]; then
    PROJECT_ROOT="$(cd "${CONDIV_PROJECT_ROOT}" && pwd)"
    SCRIPT_DIR="$PROJECT_ROOT/ConDiv_symlay"
    return 0
  fi

  if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    candidate="$(cd "${SLURM_SUBMIT_DIR}" && pwd)"
    if [ -f "$candidate/run_remote_update.sh" ] && [ -f "$candidate/layer_template.json" ]; then
      SCRIPT_DIR="$candidate"
      PROJECT_ROOT="$(cd "$candidate/.." && pwd)"
      return 0
    fi
    if [ -f "$candidate/ConDiv_symlay/run_remote_update.sh" ] && [ -f "$candidate/source.sh" ]; then
      PROJECT_ROOT="$candidate"
      SCRIPT_DIR="$candidate/ConDiv_symlay"
      return 0
    fi
  fi

  candidate="$(cd "$(dirname "$0")" && pwd)"
  if [ -f "$candidate/layer_template.json" ] && [ -f "$candidate/../source.sh" ]; then
    SCRIPT_DIR="$candidate"
    PROJECT_ROOT="$(cd "$candidate/.." && pwd)"
    return 0
  fi

  echo "ERROR: could not resolve the ConDiv_symlay workflow directory." >&2
  echo "Set CONDIV_PROJECT_ROOT or submit from the project root or ConDiv_symlay directory." >&2
  exit 1
}

resolve_paths
cd "$SCRIPT_DIR"

if [ -z "${CONDIV_ROUND_MANIFEST:-}" ]; then
  echo "ERROR: CONDIV_ROUND_MANIFEST is required."
  exit 1
fi

if [ -f /etc/profile.d/modules.sh ]; then
  source /etc/profile.d/modules.sh
fi

if command -v module >/dev/null 2>&1; then
  module load python/3.11.9
  module load cmake
  module load openmpi
fi

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OMP_PROC_BIND="${OMP_PROC_BIND:-close}"
export OMP_PLACES="${OMP_PLACES:-cores}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"
export PYTHONUNBUFFERED=1

VENV_ACTIVATE="$SCRIPT_DIR/venv/bin/activate"
if [ ! -f "$VENV_ACTIVATE" ]; then
  VENV_ACTIVATE="$PROJECT_ROOT/.venv/bin/activate"
fi

source "$VENV_ACTIVATE"
export PYTHONPATH="${PYTHONPATH:-}"
source "$PROJECT_ROOT/source.sh"

export CONDIV_PROJECT_ROOT="$PROJECT_ROOT"
export PYTHONPATH="$SCRIPT_DIR:$PROJECT_ROOT/py:$PROJECT_ROOT/obj:${PYTHONPATH:-}"

echo "Running separate-job FF update"
echo "  round manifest: ${CONDIV_ROUND_MANIFEST}"
echo "  slurm job id: ${SLURM_JOB_ID:-N/A}"

python3 -u "$SCRIPT_DIR/slurm_round.py" finalize-round --round-manifest "$CONDIV_ROUND_MANIFEST"
