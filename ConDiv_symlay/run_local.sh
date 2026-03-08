#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

source "$PROJECT_ROOT/.venv/bin/activate"
source "$PROJECT_ROOT/source.sh"

export CONDIV_PROJECT_ROOT="$PROJECT_ROOT"
export PYTHONPATH="$SCRIPT_DIR:$PROJECT_ROOT/py:$PROJECT_ROOT/obj:${PYTHONPATH:-}"

PROFILE="${PROFILE:-dimer3}"
BASE_DIR="${BASE_DIR:-$SCRIPT_DIR/test_${PROFILE}}"
RUN_STEPS="${1:-${RUN_STEPS:-20}}"
WORKER_LAUNCH="${WORKER_LAUNCH:-auto}"
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
echo "  checkpoint: $checkpoint"
echo "  steps: $RUN_STEPS"
echo "  layer manifest: $LAYER_MANIFEST_JSON"

python3 -u "$SCRIPT_DIR/ConDiv_mem.py" restart "$checkpoint" "$RUN_STEPS" 2>&1 | tee -a "$BASE_DIR/run.output"
