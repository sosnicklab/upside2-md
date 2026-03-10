#!/bin/bash
set -euo pipefail

resolve_paths() {
  local candidate

  if [ -n "${CONDIV_PROJECT_ROOT:-}" ] && [ -f "${CONDIV_PROJECT_ROOT}/source.sh" ]; then
    PROJECT_ROOT="$(cd "${CONDIV_PROJECT_ROOT}" && pwd)"
    SCRIPT_DIR="$PROJECT_ROOT/ConDiv"
    return 0
  fi

  candidate="$(cd "$(dirname "$0")" && pwd)"
  if [ -f "$candidate/ConDiv_mem.py" ] && [ -f "$candidate/../source.sh" ]; then
    SCRIPT_DIR="$candidate"
    PROJECT_ROOT="$(cd "$candidate/.." && pwd)"
    return 0
  fi

  echo "ERROR: could not resolve the ConDiv workflow directory." >&2
  echo "Run this script from ConDiv or set CONDIV_PROJECT_ROOT." >&2
  exit 1
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
BASE_DIR="$SCRIPT_DIR/validate_${PROFILE}"
INIT_FORCEFIELD_DIR="${INIT_FORCEFIELD_DIR:-$PROJECT_ROOT/parameters/ff_2.1}"
ROUNDS="${ROUNDS:-3}"
STEPS_PER_ROUND="${STEPS_PER_ROUND:-2}"
WORKER_LAUNCH="${WORKER_LAUNCH:-auto}"
FD_SAMPLES="${FD_SAMPLES:-24}"
FD_EPS="${FD_EPS:-1e-3}"
FD_REL_MEDIAN_THRESHOLD="${FD_REL_MEDIAN_THRESHOLD:-1.2e-1}"
FD_REL_MAX_THRESHOLD="${FD_REL_MAX_THRESHOLD:-6e-1}"

export CONDIV_WORKER_LAUNCH="$WORKER_LAUNCH"
export CONDIV_FF_DIR="$INIT_FORCEFIELD_DIR"

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
  python3 "$SCRIPT_DIR/ConDiv_mem.py" initialize "$INIT_FORCEFIELD_DIR" "$PROTEIN_DIR" "$PDB_LIST" "$BASE_DIR"
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

done

echo "All validation rounds completed successfully"
