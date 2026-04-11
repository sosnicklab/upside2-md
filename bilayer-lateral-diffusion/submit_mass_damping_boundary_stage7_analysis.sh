#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"

export BILAYER_DIFF_RUN_PRESET="${BILAYER_DIFF_RUN_PRESET:-mass_damping_boundary}"
export BILAYER_DIFF_BASE_DIR="${BILAYER_DIFF_BASE_DIR:-$SCRIPT_DIR/runs/$BILAYER_DIFF_RUN_PRESET}"

exec "$SCRIPT_DIR/submit_stage7_analysis.sh" "$@"
