#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
set -euo pipefail

export HYBRID_ACTIVATION_STAGE="${HYBRID_ACTIVATION_STAGE:-__hybrid_disabled__}"
export PROD_70_PREPROD_PROTEIN_MODE="${PROD_70_PREPROD_PROTEIN_MODE:-free}"

PDB_ID_VALUE="${PDB_ID:-}"
for arg in "$@"; do
    case "$arg" in
        PDB_ID=*)
            PDB_ID_VALUE="${arg#*=}"
            ;;
    esac
done

if [ -z "$PDB_ID_VALUE" ]; then
    PDB_ID_VALUE="1rkl"
fi

export RUNTIME_PDB_ID="${RUNTIME_PDB_ID:-${PDB_ID_VALUE}_nonrigid_dry}"
export RUN_DIR="${RUN_DIR:-outputs/martini_test_${PDB_ID_VALUE}_nonrigid_dry}"

exec "${SCRIPT_DIR}/run_sim_1rkl_rigid_dry.sh" "$@"
