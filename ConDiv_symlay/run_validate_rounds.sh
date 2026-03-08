#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

source "$PROJECT_ROOT/.venv/bin/activate"
source "$PROJECT_ROOT/source.sh"

export CONDIV_PROJECT_ROOT="$PROJECT_ROOT"
export PYTHONPATH="$SCRIPT_DIR:$PROJECT_ROOT/py:$PROJECT_ROOT/obj:${PYTHONPATH:-}"

PROFILE="${PROFILE:-dimer3}"
BASE_DIR="${BASE_DIR:-$SCRIPT_DIR/validate_${PROFILE}}"
INIT_FORCEFIELD_DIR="${INIT_FORCEFIELD_DIR:-$PROJECT_ROOT/parameters/ff_2.1}"
ROUNDS="${ROUNDS:-3}"
STEPS_PER_ROUND="${STEPS_PER_ROUND:-2}"
WORKER_LAUNCH="${WORKER_LAUNCH:-auto}"
FD_SAMPLES="${FD_SAMPLES:-24}"
FD_EPS="${FD_EPS:-1e-3}"
FD_REL_MEDIAN_THRESHOLD="${FD_REL_MEDIAN_THRESHOLD:-1.2e-1}"
FD_REL_MAX_THRESHOLD="${FD_REL_MAX_THRESHOLD:-6e-1}"
LAYER_MANIFEST_JSON="${LAYER_MANIFEST_JSON:-$BASE_DIR/layer_manifest.json}"

export CONDIV_WORKER_LAUNCH="$WORKER_LAUNCH"
export CONDIV_FF_DIR="$INIT_FORCEFIELD_DIR"
export CONDIV_SYMLAY_LAYER_MANIFEST="$LAYER_MANIFEST_JSON"

case "$PROFILE" in
  dimer3)
    PROTEIN_DIR="${PROTEIN_DIR:-$SCRIPT_DIR/upside_input2}"
    PDB_LIST="${PDB_LIST:-$SCRIPT_DIR/pdb_list2}"
    ;;
  legacy60)
    PROTEIN_DIR="${PROTEIN_DIR:-$SCRIPT_DIR/upside_input}"
    PDB_LIST="${PDB_LIST:-$SCRIPT_DIR/pdb_list}"
    ;;
  *)
    echo "ERROR: PROFILE must be dimer3 or legacy60"
    exit 1
    ;;
esac

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

mkdir -p "$BASE_DIR"

if [ ! -f "$BASE_DIR/initial_checkpoint.pkl" ]; then
  echo "Initializing base checkpoint from current force field"
  "$SCRIPT_DIR/run_init.sh"
fi

echo "Running multi-round validation"
echo "  profile: $PROFILE"
echo "  base dir: $BASE_DIR"
echo "  rounds: $ROUNDS"
echo "  steps/round: $STEPS_PER_ROUND"
echo "  fd samples: $FD_SAMPLES"
echo "  fd eps: $FD_EPS"

for round in $(seq 1 "$ROUNDS"); do
  echo "=== ROUND $round / $ROUNDS ==="
  checkpoint="$(find_latest_checkpoint "$BASE_DIR")"
  python3 -u "$SCRIPT_DIR/ConDiv_mem.py" restart "$checkpoint" "$STEPS_PER_ROUND" 2>&1 | tee -a "$BASE_DIR/round_${round}.log"

  checkpoint="$(find_latest_checkpoint "$BASE_DIR")"
  report="$BASE_DIR/gradient_round_${round}.json"
  python3 "$SCRIPT_DIR/check_membrane_gradient.py" \
    --checkpoint "$checkpoint" \
    --report "$report" \
    --fd-samples "$FD_SAMPLES" \
    --fd-eps "$FD_EPS" \
    --rel-median-threshold "$FD_REL_MEDIAN_THRESHOLD" \
    --rel-max-threshold "$FD_REL_MAX_THRESHOLD"

  python3 "$SCRIPT_DIR/validate_symlay_constraints.py" \
    --checkpoint "$checkpoint" \
    --layer-manifest "$LAYER_MANIFEST_JSON" \
    --report "$BASE_DIR/constraint_round_${round}.json"

done

echo "All validation rounds completed successfully"
