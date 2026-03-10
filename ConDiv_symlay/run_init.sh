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

has_existing_training_state() {
  local base_dir="$1"
  [ -f "$base_dir/initial_checkpoint.pkl" ] && return 0
  [ -f "$base_dir/training_progress.jsonl" ] && return 0
  [ -f "$base_dir/training_status.json" ] && return 0
  [ -f "$base_dir/layer_manifest.json" ] && return 0
  compgen -G "$base_dir/epoch_*_minibatch_*" > /dev/null
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
BASE_DIR="${BASE_DIR:-$SCRIPT_DIR/test_${PROFILE}}"
RUN_DIR_RECORD="$SCRIPT_DIR/.condiv_current_run_dir"
INIT_FORCEFIELD_DIR="${INIT_FORCEFIELD_DIR:-$PROJECT_ROOT/parameters/ff_2.1}"
WORKER_LAUNCH="${WORKER_LAUNCH:-auto}"
LAYER_TEMPLATE="${LAYER_TEMPLATE:-$SCRIPT_DIR/layer_template.json}"
BILAYER_PDB="${BILAYER_PDB:-$PROJECT_ROOT/example/16.MARTINI/pdb/bilayer.MARTINI.pdb}"
LIPID_ITP="${LIPID_ITP:-$PROJECT_ROOT/example/16.MARTINI/ff_dry/dry_martini_v2.1_lipids.itp}"

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

if [ -d "$BASE_DIR" ] && has_existing_training_state "$BASE_DIR"; then
  echo "ERROR: existing training run detected at $BASE_DIR"
  echo "Resume it with: sbatch run_remote.sh"
  echo "Or remove/change BASE_DIR before rerunning ./run_init.sh."
  exit 1
fi

mkdir -p "$BASE_DIR"
BASE_DIR="$(cd "$BASE_DIR" && pwd)"
LAYER_MANIFEST_JSON="${LAYER_MANIFEST_JSON:-$BASE_DIR/layer_manifest.json}"
LAYER_MANIFEST_CSV="${LAYER_MANIFEST_CSV:-$BASE_DIR/layer_manifest.csv}"
LAYER_MANIFEST_PNG="${LAYER_MANIFEST_PNG:-$BASE_DIR/layer_manifest.png}"
export CONDIV_SYMLAY_LAYER_MANIFEST="$LAYER_MANIFEST_JSON"

echo "Building ConDiv_symlay layer manifest"
echo "  layer template: $LAYER_TEMPLATE"
echo "  bilayer pdb: $BILAYER_PDB"
echo "  lipid itp: $LIPID_ITP"
echo "  python venv: $VENV_ACTIVATE"

python3 "$SCRIPT_DIR/build_layer_manifest.py" \
  --template "$LAYER_TEMPLATE" \
  --bilayer-pdb "$BILAYER_PDB" \
  --lipid-itp "$LIPID_ITP" \
  --output-json "$LAYER_MANIFEST_JSON" \
  --output-csv "$LAYER_MANIFEST_CSV" \
  --output-png "$LAYER_MANIFEST_PNG"

echo "Initializing symmetric-layer membrane ConDiv"
echo "  profile: $PROFILE"
echo "  base dir: $BASE_DIR"
echo "  init forcefield: $INIT_FORCEFIELD_DIR"
echo "  layer manifest: $LAYER_MANIFEST_JSON"

python3 "$SCRIPT_DIR/ConDiv_mem.py" initialize "$INIT_FORCEFIELD_DIR" "$PROTEIN_DIR" "$PDB_LIST" "$BASE_DIR"

python3 "$SCRIPT_DIR/validate_symlay_constraints.py" \
  --checkpoint "$BASE_DIR/initial_checkpoint.pkl" \
  --layer-manifest "$LAYER_MANIFEST_JSON" \
  --report "$BASE_DIR/constraint_init.json"

printf '%s\n' "$BASE_DIR" > "$RUN_DIR_RECORD"
echo "Recorded current run dir: $BASE_DIR"
