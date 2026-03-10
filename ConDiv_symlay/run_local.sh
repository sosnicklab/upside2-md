#!/bin/bash
set -euo pipefail

resolve_paths() {
  local candidate

  if [ -n "${CONDIV_PROJECT_ROOT:-}" ] && [ -f "${CONDIV_PROJECT_ROOT}/source.sh" ]; then
    PROJECT_ROOT="$(cd "${CONDIV_PROJECT_ROOT}" && pwd)"
    SCRIPT_DIR="$PROJECT_ROOT/ConDiv_symlay"
    return 0
  fi

  candidate="$(cd "$(dirname "$0")" && pwd)"
  if [ -f "$candidate/layer_template.json" ] && [ -f "$candidate/../source.sh" ]; then
    SCRIPT_DIR="$candidate"
    PROJECT_ROOT="$(cd "$candidate/.." && pwd)"
    return 0
  fi

  echo "ERROR: could not resolve the ConDiv_symlay workflow directory." >&2
  echo "Run this script from ConDiv_symlay or set CONDIV_PROJECT_ROOT." >&2
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
RUN_STEPS="${1:-${RUN_STEPS:-20}}"
WORKER_LAUNCH="${WORKER_LAUNCH:-auto}"

if [ -f "$RUN_DIR_RECORD" ]; then
  BASE_DIR="$(tr -d '\r' < "$RUN_DIR_RECORD")"
elif DISCOVERED_BASE_DIR="$(discover_current_checkout_run_dir)"; then
  BASE_DIR="$DISCOVERED_BASE_DIR"
else
  BASE_DIR="$SCRIPT_DIR/test_${PROFILE}"
fi
BASE_DIR="$(cd "$BASE_DIR" 2>/dev/null && pwd || printf '%s' "$BASE_DIR")"

LAYER_MANIFEST_JSON="${LAYER_MANIFEST_JSON:-$BASE_DIR/layer_manifest.json}"

export CONDIV_WORKER_LAUNCH="$WORKER_LAUNCH"
export CONDIV_SYMLAY_LAYER_MANIFEST="$LAYER_MANIFEST_JSON"

find_latest_checkpoint() {
  local base_dir="$1"
  local ckpt="$base_dir/initial_checkpoint.pkl"
  local all_dirs
  all_dirs=$(ls -d "$base_dir"/epoch_*_minibatch_* 2>/dev/null | sort -r || true)
  for d in $all_dirs; do
    if [ -f "$d/checkpoint.pkl" ]; then
      ckpt="$d/checkpoint.pkl"
      break
    fi
  done
  echo "$ckpt"
}

checkpoint="$(find_latest_checkpoint "$BASE_DIR")"
if [ ! -f "$checkpoint" ]; then
  echo "ERROR: checkpoint not found at $checkpoint"
  echo "Run ./run_init.sh first."
  exit 1
fi

echo "Running local membrane ConDiv restart"
echo "  base dir: $BASE_DIR"
if [ -f "$RUN_DIR_RECORD" ]; then
  echo "  run dir record: $RUN_DIR_RECORD"
elif [ -n "${DISCOVERED_BASE_DIR:-}" ]; then
  echo "  discovered current-checkout run dir: $DISCOVERED_BASE_DIR"
fi
echo "  checkpoint: $checkpoint"
echo "  steps: $RUN_STEPS"
echo "  layer manifest: $LAYER_MANIFEST_JSON"

python3 -u "$SCRIPT_DIR/ConDiv_mem.py" restart "$checkpoint" "$RUN_STEPS" 2>&1 | tee -a "$BASE_DIR/run.output"
