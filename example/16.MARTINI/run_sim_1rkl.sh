#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

export PDB_ID="${PDB_ID:-1rkl}"
export LIPID_RESOLUTION="${LIPID_RESOLUTION:-coarse}"
export RUNTIME_PDB_ID="${RUNTIME_PDB_ID:-1rkl_hybrid}"
export RUN_DIR="${RUN_DIR:-outputs/martini_1rkl_hybrid}"
export PROTEIN_AA_PDB="${PROTEIN_AA_PDB:-pdb/1rkl.pdb}"

exec "${SCRIPT_DIR}/run_sim_hybrid.sh" "$@"
