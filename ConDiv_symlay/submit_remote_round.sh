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
    if [ -f "$candidate/submit_remote_round.sh" ] && [ -f "$candidate/layer_template.json" ]; then
      SCRIPT_DIR="$candidate"
      PROJECT_ROOT="$(cd "$candidate/.." && pwd)"
      return 0
    fi
    if [ -f "$candidate/ConDiv_symlay/submit_remote_round.sh" ] && [ -f "$candidate/source.sh" ]; then
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
  echo "Set CONDIV_PROJECT_ROOT or run from the project root or ConDiv_symlay directory." >&2
  exit 1
}

discover_current_checkout_run_dir() {
  local candidates=()
  local ckpt
  for ckpt in "$SCRIPT_DIR"/*/initial_checkpoint.pkl; do
    [ -f "$ckpt" ] || continue
    candidates+=("$(dirname "$ckpt")")
  done
  if [ "${#candidates[@]}" -eq 0 ]; then
    return 1
  fi
  if [ "${#candidates[@]}" -eq 1 ]; then
    printf '%s' "${candidates[0]}"
    return 0
  fi
  ls -td "${candidates[@]}" | head -n 1
}

resolve_paths
cd "$SCRIPT_DIR"

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

PROFILE="${PROFILE:-dimer3}"
RUN_DIR_RECORD="$SCRIPT_DIR/.condiv_current_run_dir"
if [ "$#" -gt 1 ]; then
  echo "ERROR: submit_remote_round.sh accepts at most one optional RUN_STEPS argument."
  exit 1
fi
RUN_STEPS="${RUN_STEPS:-0}"
if [ "$#" -eq 1 ]; then
  RUN_STEPS="$1"
fi

if [ -f "$RUN_DIR_RECORD" ]; then
  RECORDED_BASE_DIR="$(tr -d '\r' < "$RUN_DIR_RECORD")"
else
  RECORDED_BASE_DIR=""
fi
if [ -n "$RECORDED_BASE_DIR" ]; then
  BASE_DIR="$RECORDED_BASE_DIR"
elif DISCOVERED_BASE_DIR="$(discover_current_checkout_run_dir)"; then
  BASE_DIR="$DISCOVERED_BASE_DIR"
else
  BASE_DIR="$SCRIPT_DIR/test_${PROFILE}"
fi
BASE_DIR="$(cd "$BASE_DIR" 2>/dev/null && pwd || printf '%s' "$BASE_DIR")"
LAYER_MANIFEST_JSON="${LAYER_MANIFEST_JSON:-$BASE_DIR/layer_manifest.json}"

PROGRESS_LOG_JSONL="${PROGRESS_LOG_JSONL:-$BASE_DIR/training_progress.jsonl}"
STATUS_JSON="${STATUS_JSON:-$BASE_DIR/training_status.json}"
CONDIV_CONVERGENCE_GRAD_NORM="${CONDIV_CONVERGENCE_GRAD_NORM:-1.0}"
CONDIV_CONVERGENCE_UPDATE_NORM="${CONDIV_CONVERGENCE_UPDATE_NORM:-0.05}"
CONDIV_CONVERGENCE_PATIENCE="${CONDIV_CONVERGENCE_PATIENCE:-3}"
CONDIV_AUTO_RESUBMIT="${CONDIV_AUTO_RESUBMIT:-1}"
CONDIV_AUTO_MAX_SUBMISSIONS="${CONDIV_AUTO_MAX_SUBMISSIONS:-0}"
CONDIV_RESUBMIT_COUNT="${CONDIV_RESUBMIT_COUNT:-0}"

export PROFILE BASE_DIR
export CONDIV_SYMLAY_LAYER_MANIFEST="$LAYER_MANIFEST_JSON"
export PROGRESS_LOG_JSONL STATUS_JSON
export CONDIV_CONVERGENCE_GRAD_NORM CONDIV_CONVERGENCE_UPDATE_NORM CONDIV_CONVERGENCE_PATIENCE
export CONDIV_AUTO_RESUBMIT CONDIV_AUTO_MAX_SUBMISSIONS CONDIV_RESUBMIT_COUNT
export RUN_STEPS

echo "Submitting separate-job Slurm round"
echo "  base dir: $BASE_DIR"
echo "  layer manifest: $LAYER_MANIFEST_JSON"
echo "  run steps remaining: $RUN_STEPS"
echo "  progress log: $PROGRESS_LOG_JSONL"
echo "  status json: $STATUS_JSON"
echo "  convergence grad threshold: ${CONDIV_CONVERGENCE_GRAD_NORM}"
echo "  convergence update threshold: ${CONDIV_CONVERGENCE_UPDATE_NORM}"
echo "  convergence patience: ${CONDIV_CONVERGENCE_PATIENCE}"
echo "  auto resubmit: ${CONDIV_AUTO_RESUBMIT}"
echo "  auto max submissions: ${CONDIV_AUTO_MAX_SUBMISSIONS}"

python3 -u "$SCRIPT_DIR/slurm_round.py" submit-round \
  --base-dir "$BASE_DIR" \
  --run-steps "$RUN_STEPS" \
  --resubmit-count "$CONDIV_RESUBMIT_COUNT"
