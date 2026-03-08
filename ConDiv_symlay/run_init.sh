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
LAYER_TEMPLATE="${LAYER_TEMPLATE:-$SCRIPT_DIR/layer_template.json}"
BILAYER_PDB="${BILAYER_PDB:-$PROJECT_ROOT/example/16.MARTINI/pdb/bilayer.MARTINI.pdb}"
LIPID_ITP="${LIPID_ITP:-$PROJECT_ROOT/example/16.MARTINI/ff_dry/dry_martini_v2.1_lipids.itp}"
LAYER_MANIFEST_JSON="${LAYER_MANIFEST_JSON:-$BASE_DIR/layer_manifest.json}"
LAYER_MANIFEST_CSV="${LAYER_MANIFEST_CSV:-$BASE_DIR/layer_manifest.csv}"
LAYER_MANIFEST_PNG="${LAYER_MANIFEST_PNG:-$BASE_DIR/layer_manifest.png}"

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

mkdir -p "$BASE_DIR"

echo "Building ConDiv_symlay layer manifest"
echo "  layer template: $LAYER_TEMPLATE"
echo "  bilayer pdb: $BILAYER_PDB"
echo "  lipid itp: $LIPID_ITP"

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
