#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"

BASE_ROOT="${BILAYER_DIFF_PROBE_BASE_ROOT:-$SCRIPT_DIR/runs/integrator_timestep_boundary}"
INTEGRATORS_CSV="${BILAYER_DIFF_PROBE_INTEGRATORS:-v,nvtc}"
TIMESTEPS_CSV="${BILAYER_DIFF_PROBE_TIME_STEPS:-0.01,0.0075,0.005}"
DAMPING_VALUES="${BILAYER_DIFF_DAMPING_VALUES:-5.0,12.0,14.0,16.0}"
MASS_SCALES="${BILAYER_DIFF_MASS_SCALES:-0.005,0.003,0.002}"
TEMPERATURE_VALUES="${BILAYER_DIFF_TEMPERATURE_VALUES:-0.70,0.90,1.10}"
REPLICATES="${BILAYER_DIFF_REPLICATES:-3}"
NVT_MAX_FORCE="${BILAYER_DIFF_PROBE_NVT_MAX_FORCE:-0.0}"

format_tag() {
  local value="$1"
  value="${value//-/m}"
  value="${value//./p}"
  printf '%s' "$value"
}

mkdir -p "$BASE_ROOT"

IFS=',' read -r -a integrators <<< "$INTEGRATORS_CSV"
IFS=',' read -r -a timesteps <<< "$TIMESTEPS_CSV"

for integrator in "${integrators[@]}"; do
  for timestep in "${timesteps[@]}"; do
    run_dir="$BASE_ROOT/${integrator}_dt$(format_tag "$timestep")"
    max_force="0.0"
    if [ "$integrator" = "nvtc" ]; then
      max_force="$NVT_MAX_FORCE"
    fi

    echo "=== submit_integrator_timestep_boundary_round: integrator=${integrator} dt=${timestep} ==="
    echo "Run directory: $run_dir"

    BILAYER_DIFF_BASE_DIR="$run_dir" \
    BILAYER_DIFF_DAMPING_VALUES="$DAMPING_VALUES" \
    BILAYER_DIFF_MASS_SCALES="$MASS_SCALES" \
    BILAYER_DIFF_TEMPERATURE_VALUES="$TEMPERATURE_VALUES" \
    BILAYER_DIFF_REPLICATES="$REPLICATES" \
    BILAYER_DIFF_TIME_STEP="$timestep" \
    BILAYER_DIFF_INTEGRATOR="$integrator" \
    BILAYER_DIFF_MAX_FORCE="$max_force" \
      "$SCRIPT_DIR/submit_remote_round.sh" "$@"
  done
done
