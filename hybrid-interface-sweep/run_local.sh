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
    candidates+=("$submit_dir" "$submit_dir/hybrid-interface-sweep")
  fi

  if [ -n "${HYBRID_SWEEP_PROJECT_ROOT:-}" ] && [ -d "${HYBRID_SWEEP_PROJECT_ROOT}" ]; then
    project_root_hint="$(cd "${HYBRID_SWEEP_PROJECT_ROOT}" && pwd)"
    candidates+=("$project_root_hint/hybrid-interface-sweep" "$project_root_hint")
  fi

  for candidate in "${candidates[@]}"; do
    if [ -f "$candidate/workflow.py" ]; then
      SCRIPT_DIR="$candidate"
      break
    fi
  done

  if [ -z "$SCRIPT_DIR" ]; then
    echo "Unable to locate the hybrid-interface-sweep workflow directory." >&2
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

source_project_env() {
  if [ -z "${PROJECT_ROOT:-}" ] || [ ! -f "$PROJECT_ROOT/source.sh" ]; then
    return 0
  fi

  if [ -z "${MY_PYTHON:-}" ]; then
    local py_path=""
    py_path="$(command -v "$PYTHON_BIN" || true)"
    if [ -n "$py_path" ]; then
      export MY_PYTHON="$(cd "$(dirname "$py_path")/.." && pwd)"
    fi
  fi

  set +u
  source "$PROJECT_ROOT/source.sh"
  set -u
}

resolve_python_bin() {
  local requested="${HYBRID_SWEEP_PYTHON:-}"
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

append_init_arg() {
  local flag="$1"
  local value="${2:-}"
  if [ -n "$value" ]; then
    INIT_CMD+=("$flag" "$value")
  fi
}

resolve_paths
cd "$SCRIPT_DIR"

if [ -f /etc/profile.d/modules.sh ]; then
  source /etc/profile.d/modules.sh
fi

if command -v module >/dev/null 2>&1; then
  module load python/3.11.9 || true
  module load cmake || true
  module load openmpi || true
  module load "${HYBRID_SWEEP_HDF5_MODULE:-hdf5/1.14.3}" || true
fi

PYTHON_BIN=""

if [ -f "$SCRIPT_DIR/.venv/bin/activate" ]; then
  source "$SCRIPT_DIR/.venv/bin/activate"
elif [ -n "${PROJECT_ROOT:-}" ] && [ -f "$PROJECT_ROOT/.venv/bin/activate" ]; then
  source "$PROJECT_ROOT/.venv/bin/activate"
fi

resolve_python_bin
source_project_env
resolve_python_bin

export PYTHONUNBUFFERED=1
if [ -n "${PROJECT_ROOT:-}" ]; then
  export HYBRID_SWEEP_PROJECT_ROOT="$PROJECT_ROOT"
fi

BASE_DIR="${HYBRID_SWEEP_BASE_DIR:-$SCRIPT_DIR/runs/default}"
mkdir -p "$BASE_DIR"

INIT_CMD=("$PYTHON_BIN" "$SCRIPT_DIR/workflow.py" init-run --base-dir "$BASE_DIR")
append_init_arg "--pdb-id" "${HYBRID_SWEEP_PDB_ID:-}"
append_init_arg "--lj-alphas" "${HYBRID_SWEEP_LJ_ALPHAS:-}"
append_init_arg "--slater-alphas" "${HYBRID_SWEEP_SLATER_ALPHAS:-}"
append_init_arg "--replicates" "${HYBRID_SWEEP_REPLICATES:-}"
append_init_arg "--seed" "${HYBRID_SWEEP_SEED:-}"
append_init_arg "--integration-ps-per-step" "${HYBRID_SWEEP_INTEGRATION_PS_PER_STEP:-}"
append_init_arg "--target-diffusion-um2-s" "${HYBRID_SWEEP_TARGET_DIFFUSION_UM2_S:-}"

if [ "${HYBRID_SWEEP_FORCE_INIT:-0}" = "1" ] || [ ! -f "$BASE_DIR/sweep_manifest.json" ]; then
  "${INIT_CMD[@]}"
fi

cmd=("$PYTHON_BIN" "$SCRIPT_DIR/workflow.py" run-local --base-dir "$BASE_DIR")
if [ -n "${HYBRID_SWEEP_MAX_TASKS:-}" ]; then
  cmd+=(--max-tasks "$HYBRID_SWEEP_MAX_TASKS")
fi
if [ -n "${HYBRID_SWEEP_START_TASK:-}" ]; then
  cmd+=(--start-task "$HYBRID_SWEEP_START_TASK")
fi
if [ "${HYBRID_SWEEP_OVERWRITE:-0}" = "1" ]; then
  cmd+=(--overwrite)
fi
if [ "${HYBRID_SWEEP_NO_ASSEMBLE:-0}" = "1" ]; then
  cmd+=(--no-assemble)
fi

"${cmd[@]}"
