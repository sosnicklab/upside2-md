#!/bin/bash
#SBATCH --job-name=condiv-mem
#SBATCH --output=slurm-%x-%j.out
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --cpus-per-task=1
#SBATCH --mem=0

set -euo pipefail

resolve_paths() {
  local candidate

  if [ -n "${CONDIV_PROJECT_ROOT:-}" ] && [ -f "${CONDIV_PROJECT_ROOT}/source.sh" ]; then
    PROJECT_ROOT="$(cd "${CONDIV_PROJECT_ROOT}" && pwd)"
    SCRIPT_DIR="$PROJECT_ROOT/ConDiv_symlay"
    return 0
  fi

  if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then
    candidate="$(cd "${SLURM_SUBMIT_DIR}" && pwd)"
    if [ -f "$candidate/run_remote.sh" ] && [ -f "$candidate/layer_template.json" ]; then
      SCRIPT_DIR="$candidate"
      PROJECT_ROOT="$(cd "$candidate/.." && pwd)"
      return 0
    fi
    if [ -f "$candidate/ConDiv_symlay/run_remote.sh" ] && [ -f "$candidate/source.sh" ]; then
      PROJECT_ROOT="$candidate"
      SCRIPT_DIR="$candidate/ConDiv_symlay"
      return 0
    fi
  fi

  candidate="$(cd "$(dirname "$0")" && pwd)"
  if [ -f "$candidate/layer_template.json" ] && [ -f "$candidate/../source.sh" ]; then
    SCRIPT_DIR="$candidate"
    PROJECT_ROOT="$(cd "$candidate/.." && pwd)"
    return 0
  fi

  echo "ERROR: could not resolve the ConDiv_symlay workflow directory." >&2
  echo "Set CONDIV_PROJECT_ROOT or submit from the project root or ConDiv_symlay directory." >&2
  exit 1
}

resolve_paths
cd "$SCRIPT_DIR"

if [ -f /etc/profile.d/modules.sh ]; then
  source /etc/profile.d/modules.sh
fi

if command -v module >/dev/null 2>&1; then
  module load python/3.11.9
  module load cmake
  module load openmpi
fi

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${VECLIB_MAXIMUM_THREADS:-1}"
export PYTHONUNBUFFERED=1

VENV_ACTIVATE="$SCRIPT_DIR/venv/bin/activate"
if [ ! -f "$VENV_ACTIVATE" ]; then
  VENV_ACTIVATE="$PROJECT_ROOT/.venv/bin/activate"
fi

source "$VENV_ACTIVATE"
export PYTHONPATH="${PYTHONPATH:-}"
source "$PROJECT_ROOT/source.sh"

export CONDIV_PROJECT_ROOT="$PROJECT_ROOT"
export PYTHONPATH="$SCRIPT_DIR:$PROJECT_ROOT/py:$PROJECT_ROOT/obj:${PYTHONPATH:-}"

PROFILE="${PROFILE:-dimer3}"
BASE_DIR="${BASE_DIR:-$SCRIPT_DIR/test_${PROFILE}}"
if [ "$#" -gt 1 ]; then
  echo "ERROR: run_remote.sh accepts at most one optional RUN_STEPS argument."
  exit 1
fi
RUN_STEPS="${RUN_STEPS:-20}"
if [ "$#" -eq 1 ]; then
  RUN_STEPS="$1"
fi
WORKER_LAUNCH="${WORKER_LAUNCH:-auto}"
LAYER_MANIFEST_JSON="${LAYER_MANIFEST_JSON:-$BASE_DIR/layer_manifest.json}"
PROGRESS_LOG_JSONL="${PROGRESS_LOG_JSONL:-$BASE_DIR/training_progress.jsonl}"
STATUS_JSON="${STATUS_JSON:-$BASE_DIR/training_status.json}"
CONDIV_CONVERGENCE_GRAD_NORM="${CONDIV_CONVERGENCE_GRAD_NORM:-1.0}"
CONDIV_CONVERGENCE_UPDATE_NORM="${CONDIV_CONVERGENCE_UPDATE_NORM:-0.05}"
CONDIV_CONVERGENCE_PATIENCE="${CONDIV_CONVERGENCE_PATIENCE:-3}"
CONDIV_AUTO_RESUBMIT="${CONDIV_AUTO_RESUBMIT:-1}"
CONDIV_AUTO_MAX_SUBMISSIONS="${CONDIV_AUTO_MAX_SUBMISSIONS:-0}"
CONDIV_RESUBMIT_COUNT="${CONDIV_RESUBMIT_COUNT:-0}"

export PROFILE BASE_DIR WORKER_LAUNCH LAYER_MANIFEST_JSON
export PROGRESS_LOG_JSONL STATUS_JSON
export CONDIV_CONVERGENCE_GRAD_NORM CONDIV_CONVERGENCE_UPDATE_NORM CONDIV_CONVERGENCE_PATIENCE
export CONDIV_AUTO_RESUBMIT CONDIV_AUTO_MAX_SUBMISSIONS CONDIV_RESUBMIT_COUNT
export RUN_STEPS

parse_slurm_task_slots() {
  local raw="${1:-}"
  raw="${raw%%,*}"
  raw="${raw%%(*}"
  printf '%s' "$raw"
}

SLURM_TASK_SLOTS="${SLURM_NTASKS:-}"
if [ -z "$SLURM_TASK_SLOTS" ]; then
  SLURM_TASK_SLOTS="$(parse_slurm_task_slots "${SLURM_NTASKS_PER_NODE:-}")"
fi

export CONDIV_N_REPLICA="${CONDIV_N_REPLICA:-8}"
export CONDIV_OMP_THREADS="${CONDIV_OMP_THREADS:-${SLURM_CPUS_PER_TASK:-1}}"
if [ -n "$SLURM_TASK_SLOTS" ] && [ "${CONDIV_N_REPLICA}" -eq "$SLURM_TASK_SLOTS" ] && [ -z "${CONDIV_MAX_PARALLEL_WORKERS:-}" ]; then
  echo "ERROR: CONDIV_N_REPLICA=$CONDIV_N_REPLICA matches all allocated Slurm task slots." >&2
  echo "This usually indicates a stale override from the older launch model." >&2
  echo "Unset CONDIV_N_REPLICA/CONDIV_N_THREADS and resubmit, or set CONDIV_MAX_PARALLEL_WORKERS explicitly if you really want one full-node replica bundle." >&2
  exit 1
fi
if [ -n "$SLURM_TASK_SLOTS" ] && [ -z "${CONDIV_MAX_PARALLEL_WORKERS:-}" ]; then
  if [ "$SLURM_TASK_SLOTS" -lt "$CONDIV_N_REPLICA" ]; then
    echo "ERROR: allocated Slurm task slots ($SLURM_TASK_SLOTS) are smaller than CONDIV_N_REPLICA ($CONDIV_N_REPLICA)." >&2
    echo "Increase --ntasks-per-node or reduce CONDIV_N_REPLICA." >&2
    exit 1
  fi
  export CONDIV_MAX_PARALLEL_WORKERS="$((SLURM_TASK_SLOTS / CONDIV_N_REPLICA))"
fi

export CONDIV_WORKER_LAUNCH="$WORKER_LAUNCH"
export CONDIV_SYMLAY_LAYER_MANIFEST="$LAYER_MANIFEST_JSON"

RESOLVED_WORKER_LAUNCH="$CONDIV_WORKER_LAUNCH"
if [ "$RESOLVED_WORKER_LAUNCH" = "auto" ] && [ -n "${SLURM_JOB_ID:-}" ]; then
  RESOLVED_WORKER_LAUNCH="srun"
fi

find_latest_checkpoint() {
  local base_dir="$1"
  local ckpt="$base_dir/initial_checkpoint.pkl"
  local all_dirs
  all_dirs=$(ls -d "$base_dir"/epoch_*_minibatch_* 2>/dev/null | sort -r || true)
  for d in $all_dirs; do
    if [ -f "$d/checkpoint.pkl" ]; then
      ckpt="$d/checkpoint.pkl"
      break
    fi
  done
  echo "$ckpt"
}

checkpoint="$(find_latest_checkpoint "$BASE_DIR")"
if [ ! -f "$checkpoint" ]; then
  echo "ERROR: checkpoint not found at $checkpoint"
  echo "Run ./run_init.sh first."
  exit 1
fi

echo "Running Slurm membrane ConDiv restart"
echo "  base dir: $BASE_DIR"
echo "  checkpoint: $checkpoint"
echo "  steps: $RUN_STEPS"
echo "  layer manifest: $LAYER_MANIFEST_JSON"
echo "  slurm job id: ${SLURM_JOB_ID:-N/A}"
echo "  upside launch request: $CONDIV_WORKER_LAUNCH"
echo "  upside launch resolved: $RESOLVED_WORKER_LAUNCH"
echo "  python venv: $VENV_ACTIVATE"
echo "  replicas/worker: ${CONDIV_N_REPLICA}"
echo "  omp threads/upside: ${CONDIV_OMP_THREADS}"
echo "  max parallel workers: ${CONDIV_MAX_PARALLEL_WORKERS:-auto}"
echo "  allocated task slots: ${SLURM_TASK_SLOTS:-unknown}"
echo "  convergence grad threshold: ${CONDIV_CONVERGENCE_GRAD_NORM}"
echo "  convergence update threshold: ${CONDIV_CONVERGENCE_UPDATE_NORM}"
echo "  convergence patience: ${CONDIV_CONVERGENCE_PATIENCE}"
echo "  auto resubmit: ${CONDIV_AUTO_RESUBMIT}"
echo "  resubmit count: ${CONDIV_RESUBMIT_COUNT}"
echo "  progress log: ${PROGRESS_LOG_JSONL}"
echo "  status json: ${STATUS_JSON}"

python3 -u "$SCRIPT_DIR/ConDiv_mem.py" restart "$checkpoint" "$RUN_STEPS" 2>&1 | tee -a "$BASE_DIR/run.output"

checkpoint="$(find_latest_checkpoint "$BASE_DIR")"
set +e
python3 "$SCRIPT_DIR/training_control.py" \
  --base-dir "$BASE_DIR" \
  --checkpoint "$checkpoint" \
  --progress-log "$PROGRESS_LOG_JSONL" \
  --status-json "$STATUS_JSON" \
  --grad-threshold "$CONDIV_CONVERGENCE_GRAD_NORM" \
  --update-threshold "$CONDIV_CONVERGENCE_UPDATE_NORM" \
  --patience "$CONDIV_CONVERGENCE_PATIENCE" \
  --run-steps "$RUN_STEPS" \
  --slurm-job-id "${SLURM_JOB_ID:-}" \
  --resubmit-count "$CONDIV_RESUBMIT_COUNT" 2>&1 | tee -a "$BASE_DIR/run.output"
control_rc=${PIPESTATUS[0]}
set -e

if [ "$control_rc" -eq 0 ]; then
  echo "Training converged; no resubmission needed."
  exit 0
fi

if [ "$control_rc" -ne 10 ]; then
  echo "ERROR: training control failed with exit code $control_rc"
  exit "$control_rc"
fi

if [ "$RUN_STEPS" -le 0 ]; then
  echo "RUN_STEPS <= 0; skipping auto-resubmission."
  exit 0
fi

if [ "$CONDIV_AUTO_RESUBMIT" != "1" ]; then
  echo "Auto-resubmission disabled; stopping before convergence."
  exit 0
fi

if [ -z "${SLURM_JOB_ID:-}" ]; then
  echo "Not running inside Slurm; skipping auto-resubmission."
  exit 0
fi

if ! command -v sbatch >/dev/null 2>&1; then
  echo "sbatch not available; skipping auto-resubmission."
  exit 0
fi

next_resubmit_count=$((CONDIV_RESUBMIT_COUNT + 1))
if [ "$CONDIV_AUTO_MAX_SUBMISSIONS" -gt 0 ] && [ "$next_resubmit_count" -gt "$CONDIV_AUTO_MAX_SUBMISSIONS" ]; then
  echo "Reached CONDIV_AUTO_MAX_SUBMISSIONS=${CONDIV_AUTO_MAX_SUBMISSIONS}; stopping before convergence."
  exit 0
fi

submit_output="$(sbatch --parsable --dependency=afterok:${SLURM_JOB_ID} --export=ALL,CONDIV_RESUBMIT_COUNT=${next_resubmit_count},RUN_STEPS=${RUN_STEPS} "$SCRIPT_DIR/run_remote.sh")"
echo "Submitted continuation job: $submit_output"
