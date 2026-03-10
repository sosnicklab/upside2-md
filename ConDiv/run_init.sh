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
BASE_DIR="$SCRIPT_DIR/test_${PROFILE}"
RUN_DIR_RECORD="$SCRIPT_DIR/.condiv_current_run_dir"
INIT_FORCEFIELD_DIR="${INIT_FORCEFIELD_DIR:-$PROJECT_ROOT/parameters/ff_2.1}"
WORKER_LAUNCH="${WORKER_LAUNCH:-auto}"

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

mkdir -p "$BASE_DIR"
BASE_DIR="$(cd "$BASE_DIR" && pwd)"

echo "Initializing membrane ConDiv"
echo "  profile: $PROFILE"
echo "  base dir: $BASE_DIR"
echo "  init forcefield: $INIT_FORCEFIELD_DIR"

python3 "$SCRIPT_DIR/ConDiv_mem.py" initialize "$INIT_FORCEFIELD_DIR" "$PROTEIN_DIR" "$PDB_LIST" "$BASE_DIR"

printf '%s\n' "$BASE_DIR" > "$RUN_DIR_RECORD"
echo "Recorded current run dir: $BASE_DIR"
