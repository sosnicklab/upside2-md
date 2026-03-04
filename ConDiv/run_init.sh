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

echo "Initializing membrane ConDiv"
echo "  profile: $PROFILE"
echo "  base dir: $BASE_DIR"
echo "  init forcefield: $INIT_FORCEFIELD_DIR"

python3 "$SCRIPT_DIR/ConDiv_mem.py" initialize "$INIT_FORCEFIELD_DIR" "$PROTEIN_DIR" "$PDB_LIST" "$BASE_DIR"
