#!/bin/bash
set -eo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/../.." && pwd)

source "$PROJECT_ROOT/.venv/bin/activate"
source "$PROJECT_ROOT/source.sh"
set -u

PDB_ID=${PDB_ID:-1rkl}
SIM_ID=${SIM_ID:-martini_hybrid}
N_REP=${N_REP:-1}
START_FRAME=${START_FRAME:-100}
PROTEIN_AA_PDB=${PROTEIN_AA_PDB:-$SCRIPT_DIR/pdb/$PDB_ID.pdb}
HYBRID_RUN_DIR=${HYBRID_RUN_DIR:-$SCRIPT_DIR/outputs/martini_${PDB_ID}_hybrid_full}
HYBRID_SOURCE_DIR=${HYBRID_SOURCE_DIR:-$HYBRID_RUN_DIR/checkpoints}
if [ -z "${HYBRID_SOURCE_PATTERN:-}" ]; then
    HYBRID_SOURCE_PATTERN="$PDB_ID.stage_7.{replica}.up"
fi
HDX_WORK_DIR=${HDX_WORK_DIR:-$HYBRID_RUN_DIR/hdx_analysis}
ANALYSIS_SCRIPT=${ANALYSIS_SCRIPT:-$PROJECT_ROOT/example/00.AnalysisScripts/analysis.sh}

WATER_ACCESSIBILITY_DIR=${WATER_ACCESSIBILITY_DIR:-}
WATER_ACCESSIBILITY_PATTERN=${WATER_ACCESSIBILITY_PATTERN:-}
HDX_T_UP=${HDX_T_UP:-}

if [ ! -f "$PROTEIN_AA_PDB" ]; then
    echo "Missing protein PDB: $PROTEIN_AA_PDB" >&2
    exit 1
fi
if [ ! -d "$HYBRID_SOURCE_DIR" ]; then
    echo "Missing hybrid trajectory directory: $HYBRID_SOURCE_DIR" >&2
    exit 1
fi

replica_token='{replica}'
first_source_name=${HYBRID_SOURCE_PATTERN//$replica_token/0}
first_source_path=$HYBRID_SOURCE_DIR/$first_source_name
if [ ! -f "$first_source_path" ]; then
    echo "Missing first hybrid trajectory: $first_source_path" >&2
    exit 1
fi
if [ -z "$HDX_T_UP" ]; then
    HDX_T_UP=$(python -c \
        'import sys, tables as tb; f = tb.open_file(sys.argv[1]); print(float(f.root.output.temperature[0, 0])); f.close()' \
        "$first_source_path")
fi

mkdir -p "$HDX_WORK_DIR/inputs" "$HDX_WORK_DIR/outputs" "$HDX_WORK_DIR/pdb" "$HDX_WORK_DIR/results"

if [ ! -f "$HDX_WORK_DIR/inputs/$PDB_ID.fasta" ]; then
    python "$UPSIDE_HOME/py/PDB_to_initial_structure.py" \
        "$PROTEIN_AA_PDB" "$HDX_WORK_DIR/inputs/$PDB_ID" --record-chain-breaks
fi
cp "$PROTEIN_AA_PDB" "$HDX_WORK_DIR/pdb/$PDB_ID.pdb"

exec env \
    runner=local \
    work_dir="$HDX_WORK_DIR" \
    pdb_id="$PDB_ID" \
    sim_id="$SIM_ID" \
    n_rep="$N_REP" \
    start_frame="$START_FRAME" \
    trajectory_mode=martini_hybrid \
    hybrid_source_dir="$HYBRID_SOURCE_DIR" \
    hybrid_source_pattern="$HYBRID_SOURCE_PATTERN" \
    hybrid_project_overwrite=true \
    water_accessibility_dir="$WATER_ACCESSIBILITY_DIR" \
    water_accessibility_pattern="$WATER_ACCESSIBILITY_PATTERN" \
    legacy_T_range="$HDX_T_UP" \
    T_targets="$HDX_T_UP" \
    target_T="$HDX_T_UP" \
    reference_temperature="$HDX_T_UP" \
    skip_experiment_data=true \
    hx_plot_enabled=auto \
    bash "$ANALYSIS_SCRIPT"
