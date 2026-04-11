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
  module load "${BILAYER_DIFF_HDF5_MODULE:-hdf5/1.14.3}" || true
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
  export BILAYER_DIFF_PROJECT_ROOT="$PROJECT_ROOT"
fi

BASE_DIR="${BILAYER_DIFF_BASE_DIR:-$SCRIPT_DIR/runs/default}"
mkdir -p "$BASE_DIR"

INIT_CMD=("$PYTHON_BIN" "$SCRIPT_DIR/workflow.py" init-run --base-dir "$BASE_DIR")
append_init_arg "--bilayer-pdb" "${BILAYER_DIFF_BILAYER_PDB:-}"
append_init_arg "--damping-values" "${BILAYER_DIFF_DAMPING_VALUES:-}"
append_init_arg "--mass-scales" "${BILAYER_DIFF_MASS_SCALES:-}"
append_init_arg "--temperature-values" "${BILAYER_DIFF_TEMPERATURE_VALUES:-}"
append_init_arg "--replicates" "${BILAYER_DIFF_REPLICATES:-}"
append_init_arg "--seed" "${BILAYER_DIFF_SEED:-}"
append_init_arg "--xy-scale" "${BILAYER_DIFF_XY_SCALE:-}"
append_init_arg "--box-padding-z" "${BILAYER_DIFF_BOX_PADDING_Z:-}"
append_init_arg "--salt-molar" "${BILAYER_DIFF_SALT_MOLAR:-}"
append_init_arg "--ion-cutoff" "${BILAYER_DIFF_ION_CUTOFF:-}"
append_init_arg "--min-60-max-iter" "${BILAYER_DIFF_MIN_60_MAX_ITER:-}"
append_init_arg "--min-61-max-iter" "${BILAYER_DIFF_MIN_61_MAX_ITER:-}"
append_init_arg "--equilibration-steps" "${BILAYER_DIFF_EQUILIBRATION_STEPS:-}"
append_init_arg "--production-steps" "${BILAYER_DIFF_PRODUCTION_STEPS:-}"
append_init_arg "--time-step" "${BILAYER_DIFF_TIME_STEP:-}"
append_init_arg "--integrator" "${BILAYER_DIFF_INTEGRATOR:-}"
append_init_arg "--max-force" "${BILAYER_DIFF_MAX_FORCE:-}"
append_init_arg "--equilibration-frame-steps" "${BILAYER_DIFF_EQUILIBRATION_FRAME_STEPS:-}"
append_init_arg "--production-frame-steps" "${BILAYER_DIFF_PRODUCTION_FRAME_STEPS:-}"

if [ "${BILAYER_DIFF_FORCE_INIT:-0}" = "1" ] || [ ! -f "$BASE_DIR/diffusion_manifest.json" ]; then
  "${INIT_CMD[@]}"
fi

cmd=("$PYTHON_BIN" "$SCRIPT_DIR/workflow.py" run-local --base-dir "$BASE_DIR")
if [ -n "${BILAYER_DIFF_MAX_TASKS:-}" ]; then
  cmd+=(--max-tasks "$BILAYER_DIFF_MAX_TASKS")
fi
if [ -n "${BILAYER_DIFF_START_TASK:-}" ]; then
  cmd+=(--start-task "$BILAYER_DIFF_START_TASK")
fi
if [ "${BILAYER_DIFF_OVERWRITE:-0}" = "1" ]; then
  cmd+=(--overwrite)
fi
if [ "${BILAYER_DIFF_NO_ASSEMBLE:-0}" = "1" ]; then
  cmd+=(--no-assemble)
fi

"${cmd[@]}"
