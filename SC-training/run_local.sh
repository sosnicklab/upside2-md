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
    candidates+=("$submit_dir" "$submit_dir/SC-training")
  fi

  if [ -n "${SC_TRAIN_PROJECT_ROOT:-}" ] && [ -d "${SC_TRAIN_PROJECT_ROOT}" ]; then
    project_root_hint="$(cd "${SC_TRAIN_PROJECT_ROOT}" && pwd)"
    candidates+=("$project_root_hint/SC-training" "$project_root_hint")
  fi

  for candidate in "${candidates[@]}"; do
    if [ -f "$candidate/workflow.py" ]; then
      SCRIPT_DIR="$candidate"
      break
    fi
  done

  if [ -z "$SCRIPT_DIR" ]; then
    echo "Unable to locate the SC-training workflow directory." >&2
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

resolve_paths
cd "$SCRIPT_DIR"

PYTHON_BIN="${SC_TRAIN_PYTHON:-python3}"

if [ -f "$SCRIPT_DIR/.venv/bin/activate" ]; then
  source "$SCRIPT_DIR/.venv/bin/activate"
elif [ -n "$PROJECT_ROOT" ] && [ -f "$PROJECT_ROOT/.venv/bin/activate" ]; then
  source "$PROJECT_ROOT/.venv/bin/activate"
fi

source_project_env

export PYTHONUNBUFFERED=1
if [ -n "$PROJECT_ROOT" ]; then
  export SC_TRAIN_PROJECT_ROOT="$PROJECT_ROOT"
fi

BASE_DIR="${SC_TRAIN_BASE_DIR:-$SCRIPT_DIR/runs/default}"
mkdir -p "$BASE_DIR"

if [ ! -f "$BASE_DIR/training_manifest.json" ]; then
  "$PYTHON_BIN" "$SCRIPT_DIR/workflow.py" init-run --base-dir "$BASE_DIR"
fi

"$PYTHON_BIN" "$SCRIPT_DIR/workflow.py" run-local --base-dir "$BASE_DIR"

if [ "${RUN_BENCHMARK:-0}" = "1" ]; then
  "$PYTHON_BIN" "$SCRIPT_DIR/workflow.py" run-benchmark --base-dir "$BASE_DIR" --execute
fi
