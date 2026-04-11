#!/bin/bash
set -euo pipefail

resolve_paths() {
  local candidate=""
  local self_dir=""
  local submit_dir=""
  local project_root_hint=""
  local -a candidates=()
  self_dir="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
  SCRIPT_DIR=""
  PROJECT_ROOT=""

  candidates+=("$self_dir")

  if [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -d "${SLURM_SUBMIT_DIR}" ]; then
    submit_dir="$(cd "${SLURM_SUBMIT_DIR}" && pwd)"
    candidates+=("$submit_dir" "$submit_dir/bilayer-lateral-diffusion")
  fi

  if [ -n "${BILAYER_DIFF_PROJECT_ROOT:-}" ] && [ -d "${BILAYER_DIFF_PROJECT_ROOT}" ]; then
    project_root_hint="$(cd "${BILAYER_DIFF_PROJECT_ROOT}" && pwd)"
    candidates+=("$project_root_hint/bilayer-lateral-diffusion" "$project_root_hint")
  fi

  for candidate in "${candidates[@]}"; do
    if [ -f "$candidate/workflow.py" ]; then
      SCRIPT_DIR="$candidate"
      break
    fi
  done

  if [ -z "$SCRIPT_DIR" ]; then
    echo "Unable to locate the bilayer-lateral-diffusion workflow directory." >&2
    return 1
  fi

  if [ -f "$SCRIPT_DIR/../source.sh" ]; then
    PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
  elif [ -n "$project_root_hint" ] && [ -f "$project_root_hint/source.sh" ]; then
    PROJECT_ROOT="$project_root_hint"
  elif [ -n "$submit_dir" ] && [ -f "$submit_dir/source.sh" ]; then
    PROJECT_ROOT="$submit_dir"
  fi
}

resolve_python_bin() {
  local requested="${BILAYER_DIFF_PYTHON:-}"
  if [ -n "$requested" ]; then
    if command -v "$requested" >/dev/null 2>&1; then
      PYTHON_BIN="$(command -v "$requested")"
    else
      PYTHON_BIN="$requested"
    fi
    return
  fi

  if [ -n "${VIRTUAL_ENV:-}" ] && [ -x "${VIRTUAL_ENV}/bin/python3" ]; then
    PYTHON_BIN="${VIRTUAL_ENV}/bin/python3"
    return
  fi

  PYTHON_BIN="$(command -v python3)"
}

resolve_paths
cd "$SCRIPT_DIR"

if [ -f /etc/profile.d/modules.sh ]; then
  source /etc/profile.d/modules.sh
fi

if command -v module >/dev/null 2>&1; then
  module load python/3.11.9 || true
fi

PYTHON_BIN=""

if [ -f "$SCRIPT_DIR/.venv/bin/activate" ]; then
  source "$SCRIPT_DIR/.venv/bin/activate"
elif [ -n "${PROJECT_ROOT:-}" ] && [ -f "$PROJECT_ROOT/.venv/bin/activate" ]; then
  source "$PROJECT_ROOT/.venv/bin/activate"
fi

resolve_python_bin

BASE_DIR="${BILAYER_DIFF_BASE_DIR:-$SCRIPT_DIR/runs/default}"
cmd=("$PYTHON_BIN" "$SCRIPT_DIR/workflow.py" pack-stage7-results --base-dir "$BASE_DIR")
if [ -n "${BILAYER_DIFF_PACK_OUTPUT:-}" ]; then
  cmd+=(--output "$BILAYER_DIFF_PACK_OUTPUT")
fi
if [ "${BILAYER_DIFF_PACK_METADATA_ONLY:-0}" = "1" ]; then
  cmd+=(--metadata-only)
fi

"${cmd[@]}"
