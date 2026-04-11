#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"

BASE_ROOT="${BILAYER_DIFF_PROBE_BASE_ROOT:-$SCRIPT_DIR/runs/integrator_timestep_boundary}"
INTEGRATORS_CSV="${BILAYER_DIFF_PROBE_INTEGRATORS:-v,nvtc}"
TIMESTEPS_CSV="${BILAYER_DIFF_PROBE_TIME_STEPS:-0.01,0.0075,0.005}"

format_tag() {
  local value="$1"
  value="${value//-/m}"
  value="${value//./p}"
  printf '%s' "$value"
}

has_stage7_output() {
  local run_dir="$1"
  if [ ! -d "$run_dir/tasks" ]; then
    return 1
  fi
  find "$run_dir/tasks" -name 'bilayer_diffusion.stage_7.0.up' -print -quit | grep -q .
}

IFS=',' read -r -a integrators <<< "$INTEGRATORS_CSV"
IFS=',' read -r -a timesteps <<< "$TIMESTEPS_CSV"

for integrator in "${integrators[@]}"; do
  for timestep in "${timesteps[@]}"; do
    run_dir="$BASE_ROOT/${integrator}_dt$(format_tag "$timestep")"
    if ! has_stage7_output "$run_dir"; then
      echo "Skipping stage-7 analysis for $run_dir: no stage_7.0 outputs found"
      continue
    fi

    echo "=== submit_integrator_timestep_boundary_stage7_analysis: integrator=${integrator} dt=${timestep} ==="
    echo "Run directory: $run_dir"

    BILAYER_DIFF_BASE_DIR="$run_dir" \
      "$SCRIPT_DIR/submit_stage7_analysis.sh" "$@"
  done
done
